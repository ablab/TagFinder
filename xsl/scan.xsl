<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>


    <xsl:template match="scan">
        <html>
            <title>Spectrum graph</title>
            <body>
                <xsl:apply-templates select="table"/>
            </body>
        </html>
    </xsl:template>

    <xsl:template match="table">
        <h4>Component <xsl:value-of select="position()"/></h4>
        <table>
            <xsl:apply-templates select="row"/>
        </table>
    </xsl:template>

    <xsl:template match="row">
        <tr>
            <xsl:apply-templates select="cell"/>
        </tr>
    </xsl:template>

    <xsl:template match="cell">
        <td>
            <xsl:apply-templates/>
        </td>
    </xsl:template>

    <xsl:template match="peak">
        <a href="#" title="{value}">
            <xsl:choose>
                <xsl:when test="type = 'B'">&gt;</xsl:when>
                <xsl:otherwise>&lt;</xsl:otherwise>
            </xsl:choose>
        </a>
    </xsl:template>

    <xsl:template match="acid">
        <xsl:value-of select="name"/>
    </xsl:template>

</xsl:stylesheet>