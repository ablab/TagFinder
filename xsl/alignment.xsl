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
                <div>Precursor mass <xsl:value-of select="scan/precursor-mass"/></div>
                <div>Protein mass <xsl:value-of select="protein/protein-mass"/></div>

                <script>
                    var width = 50;
                    function drawScale(ctx) {
                        ctx.strokeStyle = "red";
                        ctx.fillStyle = "red";
                        for (var i = -15; 35 >=i; i++) {
                            ctx.beginPath();
                            x = 250 + i * width;
                            ctx.moveTo(x, 500);
                            ctx.lineTo(x, 450);
                            ctx.stroke();
                            ctx.fillText(i, x +3, 470);
                        }
                    }
                    var peaks = [<xsl:value-of select="peaks"/>];                    
                    
                </script>
                
                <span style="font-family:monospace">
                <xsl:value-of select="protein/protein-sequence"/>

                    <br/>
                    <xsl:apply-templates select="cleavage" mode="line"/>
                    <br/>
                    Scale <input type="range" id="scale" min="0" max="3000" value="100" style="width:500;" onchange="drawPeaks()"/><span id="scaleDisplay"></span>3000%
                    <a href="#" onclick="document.getElementById('scale').value=100;">Reset</a>
                    <br/>

                <canvas id="canvas" width="1000px" height="500"/>
<script>
    var ctx = document.getElementById("canvas").getContext('2d');

    var target;
    var emass;
    var modification;

    function drawPeaks(targetParam, emassParam, modificationParam) {
        if (targetParam != null) {
            target = targetParam;
        }
        if (emassParam == null) {
            emass = [];
        }
        if (modificationParam == null) {
            modification = 0;
        }
        var scale = document.getElementById('scale').value;
        var x;
        ctx.fillStyle = "white";
        ctx.clearRect(0, 0, 100000, 100000);
        ctx.strokeStyle = "black";
        ctx.fillStyle = "black";
        ctx.fillText("monomass : " + target, 800, 50);
        ctx.fillText("scale : " + scale + "%", 825, 60);
        drawScale(ctx);

        ctx.strokeStyle = "green";
        for (var i = 0; emass.length > i; i++) {
            ctx.beginPath();
            x = 249 + (emass[i][0] + modification - target) * width;
            ctx.moveTo(x, 500);
            ctx.lineTo(x, 500 - emass[i][1] * 5);
            ctx.stroke();
        }
        
        ctx.strokeStyle = "black";


        for (var i = 0; peaks.length > i; i++) {
            var m = peaks[i][0];
            if (15 > Math.abs(m - target)) {
                ctx.beginPath();
                x = 251 + (m - target) * width;
                ctx.moveTo(x, 500);
                ctx.lineTo(x, 500 - peaks[i][1] / 1000 * scale / 100);
                ctx.stroke();
            }
        }
    }
</script>
                    <br/>
                    <xsl:apply-templates select="scan/peaks/peak" mode="deconv"/>
                <br/>
                <xsl:apply-templates select="cleavage"/>
                </span>

            </body>
        </html>
    </xsl:template>

    <xsl:template match="cleavage" mode="line">
        <xsl:apply-templates select="peptide"/>
        <br/>
    </xsl:template>

    <xsl:template match="peptide">
        <a href="#" onclick="drawPeaks({mass}, emassc{../position}_{position()},{../modification}); return false"><xsl:value-of select="acid"/></a>
        <script>
            var emassc<xsl:value-of select="../position"/>_<xsl:value-of select="position()"/> = [
            <xsl:apply-templates select="tpeak" mode="data"/>
            ]
        </script>
    </xsl:template>

    <xsl:template match="peak" mode="deconv">
        <a href="#" onclick="drawPeaks({mass}, null, null); return false"><xsl:value-of select="mass"/>&#160;<xsl:value-of select="intencity"/></a> <br/>
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
                <xsl:value-of select="substring(../../protein/protein-sequence, end-position  + 1, $pos - (end-position)  +1)"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:call-template name="ident">
                    <xsl:with-param name="pos" select="$pos"/>
                </xsl:call-template>
                <xsl:value-of select="substring(../../protein/protein-sequence, $pos  + 1, (end-position) - $pos)"/>
            </xsl:otherwise>
        </xsl:choose>

        <xsl:apply-templates select="modification"/> <a href="#" onclick="drawPeaks({expected-mass}, emass{$pos}_{position()}, {protein-modification});return false;">draw</a>
        <script>
            var emass<xsl:value-of select="$pos"/>_<xsl:value-of select="position()"/> = [
                <xsl:apply-templates select="tpeak" mode="data"/>
            ]
        </script>



        <!--<canvas id="canvas{$pos}_{position()}" width="1000px" height="500"/>-->
        <!--<script>-->
            <!--var ctx = document.getElementById("canvas<xsl:value-of select="$pos"/>_<xsl:value-of select="position()"/>").getContext('2d');-->
            <!--var x;-->
            <!--var target = ;-->
            <!--drawScale(ctx);-->
            <!--ctx.strokeStyle = "green";-->
            <!--<xsl:apply-templates select="tpeak"/>-->
            <!--ctx.strokeStyle = "black";-->
            <!--&lt;!&ndash;<xsl:apply-templates select="peak/mz-support" mode="draw"/>&ndash;&gt;-->
            <!--for (var i = 0; peaks.length > i; i++) {-->
                <!--var m = peaks[i][0];-->
                <!--if (15 > Math.abs(m - target)) {-->
                    <!--ctx.beginPath();-->
                    <!--var x = 501 + (m - target) * 30;-->
                    <!--ctx.moveTo(x, 500);-->
                    <!--ctx.lineTo(x, 500 - peaks[i][1]/1000);-->
                    <!--ctx.stroke();-->
                <!--}-->
            <!--}-->

        <!--</script>-->
        <!--<xsl:apply-templates select="peak"/>-->
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

    <xsl:template match="peak">
        <xsl:value-of select="delta"/>
        <xsl:apply-templates select="mz-support"/>
        <br/>
    </xsl:template>

    <xsl:template match="mz-support">
        &#160;<xsl:value-of select="charge"/>&#160;<xsl:value-of select="intencity"/>&#160;<xsl:value-of select="error"/> &#160;<xsl:value-of select="mass * charge"/>
    </xsl:template>

    <xsl:template match="mz-support" mode="draw">
        ctx.beginPath();
        x = 501 + (<xsl:value-of select="mass * charge"/> - <xsl:value-of select="charge"/> - target) * 30;
        ctx.moveTo(x, 500);
        ctx.lineTo(x, 500 - <xsl:value-of select="intencity"/>/1000);
        ctx.stroke();
    </xsl:template>

    <xsl:template match="tpeak">
        ctx.beginPath();
        x = 499 + (<xsl:value-of select="value"/> - target + <xsl:value-of select="../protein-modification"/>) * 30;
        ctx.moveTo(x, 500);
        ctx.lineTo(x, 500 - <xsl:value-of select="score"/> * 5);
        ctx.stroke();
    </xsl:template>

    <xsl:template match="tpeak" mode="data">
        [<xsl:value-of select="value"/>, <xsl:value-of select="score"/>],
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