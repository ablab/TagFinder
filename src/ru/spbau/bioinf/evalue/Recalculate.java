package ru.spbau.bioinf.evalue;

import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import org.apache.log4j.Logger;
import ru.spbau.bioinf.tagfinder.Configuration;
import ru.spbau.bioinf.tagfinder.Scan;

public class Recalculate {

    private static Logger log = Logger.getLogger(Recalculate.class);

    public static void main(String[] args) throws Exception {
        Configuration conf = new Configuration(args);
        EValueServer.init(args);
        Map<Integer,Scan> scans = conf.getScans();
        DbUtil.initDatabase();
        Connection con = DbUtil.getConnection();
        PreparedStatement ps = null;
        ResultSet rs = null;
        List<int[]> process = new ArrayList<int[]>();
        try {
            ps = con.prepareStatement("select scan_id, protein_id from t_history order by evalue");
            rs = ps.executeQuery();
            while (rs.next()) {
                process.add(new int[]{rs.getInt(1), rs.getInt(2)});
            }
        } catch (SQLException e) {
            log.error("Error loading data from database", e);
            throw new RuntimeException(e);
        } finally {
            DbUtil.close(con, ps, rs);
        }

        for (int[] pair : process) {
                EValueServer.getEvalue(pair[0], pair[1]);
        }
    }
}
