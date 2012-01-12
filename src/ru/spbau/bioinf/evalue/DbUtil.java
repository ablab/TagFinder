package ru.spbau.bioinf.evalue;

import java.sql.Connection;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import org.apache.log4j.Logger;
import org.postgresql.ds.PGPoolingDataSource;

public class DbUtil {

    private static Logger log = Logger.getLogger(DbUtil.class);
    
    private static final PGPoolingDataSource dataSource = new PGPoolingDataSource();

    public static void initDatabase() {
        dataSource.setDataSourceName("EValue database");
        dataSource.setServerName("127.0.0.1");
        dataSource.setDatabaseName("salmonella");
        dataSource.setUser("postgres");
        dataSource.setPassword("yasha");
        dataSource.setMaxConnections(10);
    }

    public static Connection getConnection() {
        try {
            return dataSource.getConnection();
        } catch (SQLException e) {
            log.error("Database failed", e);
            throw new RuntimeException(e);
        }
    }

    public static void close(Connection con) {
        try {
            con.close();
        } catch (SQLException e) {
            log.error("Error closing connection", e);
        }
    }

    public static void close(Connection con, Statement st) {
        close(st);
        close(con);
    }

    public static void close(Connection con, Statement st, ResultSet rs) {
        close(rs);
        close(st);
        close(con);
    }

    public static void close(Statement st) {
        if (st != null) {
            try {
                st.close();
            } catch (SQLException e) {
                log.error("Error closing statement", e);
            }
        }
    }

    public static void close(ResultSet rs) {
        if (rs != null) {
            try {
                rs.close();
            } catch (SQLException e) {
                log.error("Error closing ResultSet", e);
            }
        }
    }
}
