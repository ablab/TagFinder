package ru.spbau.bioinf.evalue;


import edu.ucsd.msalign.align.prsm.PrSM;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.List;
import java.util.Map;
import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;
import javax.xml.transform.stream.StreamSource;
import net.sf.saxon.s9api.Processor;
import net.sf.saxon.s9api.Serializer;
import net.sf.saxon.s9api.XdmNode;
import net.sf.saxon.s9api.XsltCompiler;
import net.sf.saxon.s9api.XsltTransformer;
import org.apache.log4j.Logger;
import org.eclipse.jetty.server.Request;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.AbstractHandler;
import org.jdom.Element;
import org.jdom.transform.JDOMSource;
import ru.spbau.bioinf.palign.Aligner;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.EValueAdapter;
import ru.spbau.bioinf.tagfinder.Protein;
import ru.spbau.bioinf.tagfinder.Scan;


public class EValueServer extends AbstractHandler {
    private static Map<Integer, Scan> scans;

    private static Logger log = Logger.getLogger(EValueServer.class);
    private static List<Protein> proteins;

    public static void main(String[] args) throws Exception {
        init(args);

        Server server = new Server(8080);
        server.setHandler(new EValueServer());
        server.start();
    }

    public static void init(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        EValueAdapter.init(conf);
        scans = conf.getScans();
        proteins = conf.getProteins();
        DbUtil.initDatabase();
    }

    public void handle(String s, Request request, HttpServletRequest httpServletRequest, HttpServletResponse response) throws IOException, ServletException {
        PrintWriter out = response.getWriter();
        try {
            int scanId = Integer.parseInt(httpServletRequest.getParameter("scanId"));
            int proteinId = Integer.parseInt(httpServletRequest.getParameter("proteinId"));
            String path = request.getPathInfo();
            if (path.endsWith("evalue")) {
                out.println(Double.toString(getEvalue(scanId, proteinId)));
            } else if (path.endsWith("align")) {
                Element alignment = Aligner.findAlignment(scans.get(scanId), proteins.get(proteinId)).toXml();
                Processor processor = new Processor(false);
                XdmNode source = processor.newDocumentBuilder().build(new JDOMSource(alignment));
                Serializer ser = new Serializer();
                ser.setOutputProperty(Serializer.Property.METHOD, "html");
                ser.setOutputProperty(Serializer.Property.INDENT, "yes");
                ser.setOutputWriter(out);

                XsltCompiler comp = processor.newXsltCompiler();
                XsltTransformer trans = comp.compile(new StreamSource(
                        new InputStreamReader(new FileInputStream(new File("xsl", "alignment.xsl")), "UTF-8"))).load();
                trans.setInitialContextNode(source);
                trans.setDestination(ser);
                trans.transform();

                //XmlUtil.outputter.output(alignment, out);
            }
        } catch (Throwable e) {
            log.debug("Error executing request", e);
            out.println("Error: " + e.getMessage());
        } finally {
            out.close();
        }
    }


    public static double getEvalue(int scanId, int proteinId) throws Exception {
        Connection con = DbUtil.getConnection();
        PreparedStatement ps = null;
        ResultSet rs = null;
        try {
            ps = con.prepareStatement("select evalue from t_status where scan_id = ? and protein_id = ?");
            ps.setInt(1, scanId);
            ps.setInt(2, proteinId);
            rs = ps.executeQuery();
            if (rs.next()) {
                return rs.getDouble(1);
            }
        } catch (SQLException e) {
            log.error("Error while loading E-value from database", e);
            throw new RuntimeException(e);
        } finally {
            DbUtil.close(con, ps, rs);
        }
        PrSM bestEValue = EValueAdapter.getBestEValue(scans.get(scanId), proteinId);
        double evalue = bestEValue != null ? bestEValue.getEValue() : 9E100;
        con = DbUtil.getConnection();
        try {
            ps = con.prepareStatement("insert into t_history (scan_id, protein_id, evalue, version) values(?, ?, ?, ?)");
            ps.setInt(1, scanId);
            ps.setInt(2, proteinId);
            ps.setDouble(3, evalue);
            ps.setString(4, EValueAdapter.getVersion());
            ps.execute();
            DbUtil.close(ps);
            ps = con.prepareStatement("insert into t_status (scan_id, protein_id, evalue) values(?, ?, ?)");
            ps.setInt(1, scanId);
            ps.setInt(2, proteinId);
            ps.setDouble(3, evalue);
            ps.execute();
        } catch (SQLException e) {
            log.error("Error while saving E-value to database", e);
            throw new RuntimeException(e);
        } finally {
            DbUtil.close(con, ps);
        }

        return evalue;
    }

}
