<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">
    <xsl:output encoding="UTF-8" method="html" omit-xml-declaration="yes" indent="yes"/>

    <xsl:variable name="width" select="1900"/>

    <xsl:template match="prsm">
        <html>
            <title>Scan #<xsl:value-of select="scan/scan-id"/> matched to protein #<xsl:value-of select="protein/protein-id"/></title>
            <script type="text/javascript" src="../../../../js/prsm.js"/>
            <body>
                <h3>Scan #<xsl:value-of select="scan/scan-id"/></h3>
                <div>Precursor Mass <xsl:value-of select="scan/precursor-mass"/>, <xsl:value-of select="count(scan/peaks/peak)"/> peaks. </div>
                <h3>Protein #<xsl:value-of select="protein/protein-id"/>
                    <xsl:text> </xsl:text>
                    <xsl:value-of select="protein/protein-name"/></h3>
                <script>
                    var sequence = '<xsl:value-of select="protein/protein-sequence"/>';
                </script>

                <div>
                Scale <input id="scale" value="100" onchange="repaintPrsm();"/>%   <a onclick="repaintPrsm();" href="#">Repaint</a>
                </div>

                <div id="prefix">...</div>


                <div id="canvasesdiv" style="position:relative; width:{$width}px; height:1000px">
                    <canvas id="prsm" style="z-index: 1; position:absolute; left:0px; top:0px;" width="{$width}" height="1000"/>
                    <canvas id="tooltip" style="z-index: 2; position:absolute; left:0px; top:0px;" width="{$width}" height="1000"/>
                </div>

                <textarea id="deltas" cols="50" rows="5"></textarea>

                <script>
                    var prsmCanvas = document.getElementById('prsm');
                    var tooltipCanvas = document.getElementById('tooltip');
                    var ctx = prsmCanvas.getContext('2d');
                    var tCtx = tooltipCanvas.getContext('2d');
                    var font = "10pt Arial";
                    ctx.font = font;
                    var peaks = [
                        <xsl:apply-templates select="scan/peaks/peak"/>
                    ];


                    var prefixLen = 0;
                    initPrsm();

                    document.getElementById('deltas').value = deltas;
                    var prefix = document.getElementById('prefix');
                    repaintPrsm();

                    var scaleControl =  document.getElementById('scale');
                    function doKeyDown(e) {
                        switch (e.keyCode) {
                            case 107: scaleControl.value++; repaintPrsm(); break;
                            case 109: document.getElementById('scale').value--; repaintPrsm(); break;
                            case 39: if (sequence.length > prefixLen) prefixLen++; leftRightMove(); break;
                            case 37: if (prefixLen > 0) prefixLen--; leftRightMove(); break;
                        }
                    }

                    function leftRightMove() {
                        prefix.innerHTML = sequence.substr(0, prefixLen) + "...";
                        repaintPrsm();
                    }

                    window.addEventListener('keydown', doKeyDown, true);
                </script>
            </body>
        </html>
    </xsl:template>


    <xsl:template match="peak">
        {mass: <xsl:value-of select="mass"/>,
        intencity: <xsl:value-of select="intencity"/>,
        charge: <xsl:value-of  select="charge"/>},
    </xsl:template>


</xsl:stylesheet>