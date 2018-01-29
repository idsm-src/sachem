/*
 * Copyright (C) 2015-2017 Jakub Galgonek   galgonek@uochb.cas.cz
 *
 * This program is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either version 2.1 of the License, or (at your option)
 * any later version. All we ask is that proper credit is given for our work, which includes - but is not limited to -
 * adding the above copyright notice to the beginning of your source code files, and to any copyright notice that you
 * may distribute with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
 * details.
 */
package cz.iocb.sachem.load;

import java.sql.Connection;
import java.sql.SQLException;
import java.util.Properties;
import org.postgresql.ds.PGPoolingDataSource;



public class ConnectionPool
{
    private PGPoolingDataSource connectionPool = null;


    public ConnectionPool(Properties properties)
    {
        connectionPool = new PGPoolingDataSource();
        connectionPool.setDataSourceName("PubChem Data Source");
        connectionPool.setServerName(properties.getProperty("host"));
        connectionPool.setPortNumber(Integer.parseInt(properties.getProperty("port")));
        connectionPool.setDatabaseName(properties.getProperty("database"));
        connectionPool.setUser(properties.getProperty("user"));
        connectionPool.setPassword(properties.getProperty("password"));
        connectionPool.setSocketTimeout(Integer.parseInt(properties.getProperty("socketTimeout")));
        connectionPool.setTcpKeepAlive(properties.getProperty("tcpKeepAlive").equals("true"));
        connectionPool.setCompatible(properties.getProperty("assumeMinServerVersion"));
        connectionPool.setMaxConnections(Integer.parseInt(properties.getProperty("maxConnections")));
    }


    public Connection getConnection() throws SQLException
    {
        return connectionPool.getConnection();
    }
}
