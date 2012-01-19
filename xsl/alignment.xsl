<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:variable name="width">500</xsl:variable>

    <xsl:template match="alignment">
        <html>
            <title>Alignment for scan #<xsl:value-of select="scan/scan-id"/> and protein #<xsl:value-of select="protein/protein-id"/></title>
            <body>
                <h3>Alignment for scan #<xsl:value-of select="scan/scan-id"/> and protein #<xsl:value-of select="protein/protein-id"/></h3>

                <div><xsl:value-of select="peaks-described"/> peaks described out of <xsl:value-of select="peaks-total"/> peaks total.</div>

                <span style="font-family:monospace">
                <xsl:value-of select="protein/protein-sequence"/>
                <br/>
                <xsl:apply-templates select="cleavage"/>
                </span>

            </body>
        </html>
    </xsl:template>

    <xsl:template match="cleavage">
        <xsl:value-of select="modification"/><br/>
        <xsl:apply-templates select="support">
            <xsl:with-param name="pos" select="position"/>
        </xsl:apply-templates>
    </xsl:template>

    <xsl:template match="support">
        <xsl:param name="pos"/>
        <xsl:choose>
            <xsl:when test="$pos > end-position">
                <xsl:call-template name="ident">
                    <xsl:with-param name="pos" select="end-position"/>
                </xsl:call-template>
                <xsl:call-template name="acid-mark">
                    <xsl:with-param name="pos" select="$pos - (end-position)"/>
                </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>
                <xsl:call-template name="ident">
                    <xsl:with-param name="pos" select="$pos"/>
                </xsl:call-template>
                <xsl:call-template name="acid-mark">
                    <xsl:with-param name="pos" select="end-position - $pos"/>
                </xsl:call-template>
            </xsl:otherwise>
        </xsl:choose>

        <xsl:apply-templates select="modification"/>
        <br/>

    </xsl:template>

    <xsl:template match="modification">
        <a href="#" onclick="return false;" style="text-decoration:none" title="{peak/mass} {error}">
        <xsl:choose>
            <xsl:when test="modification-type='NONE'"><font color="green">!</font></xsl:when>
            <xsl:when test="modification-type='ONE_LOSS'"><font color="red">-1</font></xsl:when>
            <xsl:when test="modification-type='ONE_GAIN'"><font color="red">+1</font></xsl:when>
            <xsl:when test="modification-type='WATER_LOSS'"><font color="blue">-W</font></xsl:when>
            <xsl:otherwise><xsl:value-of select="modification-type"/></xsl:otherwise>
        </xsl:choose>
        </a>
    </xsl:template>


    <xsl:template name="ident">
        <xsl:param name="pos"/>
        <xsl:if test="$pos > 0">
<xsl:text>&#160;</xsl:text>
            <xsl:call-template name="ident">
                <xsl:with-param name="pos" select="($pos) - 1"/>
            </xsl:call-template>
        </xsl:if>
    </xsl:template>

    <xsl:template name="acid-mark">
        <xsl:param name="pos"/>
        <xsl:if test="$pos > 0">
            <xsl:text>^</xsl:text>
            <xsl:call-template name="acid-mark">
                <xsl:with-param name="pos" select="($pos) - 1"/>
            </xsl:call-template>
        </xsl:if>
    </xsl:template>

</xsl:stylesheet>