<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:output method="html" encoding="iso-8859-1" indent="yes"
 doctype-public="-//W3C//DTD XHTML 1.0 Transitional//EN"
 doctype-system="http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd"/>

<xsl:template match="/CP2K_INPUT">
 <html xmlns="http://www.w3.org/1999/xhtml" lang="en" xml:lang="en">
  <head>
   <title>CP2K input</title>
   <style type="text/css">
    body {background-color: #eeeeee}
    big.uctt {font-family: monospace; text-transform: uppercase}
    ul.none {list-style-type: none}
    ul.disc {list-style-type: disc}
    ul.circle {list-style-type: circle}
   </style>
  </head>
  <body>
   <h1 align="center">CP2K input description</h1>
   <h2>Version information</h2>
    This HTML manual was generated automatically from a CP2K executable
    compiled on <xsl:value-of select="COMPILE_DATE"/> using the
    --xml command line option. Thus the manual describes exactly this
    version of the code. The latest CVS log file entry found was
    <xsl:value-of select="COMPILE_LASTCVS"/>.
   <h2>Input structure</h2>
    All sections that can be part of a CP2K input file are shown
    with their allowed nestings. A detailed description can be obtained by
    clicking on the links. The links in the detailed descriptions switch
    back to the corresponding index entry. In this way a toggling between
    the index and the detailed description is feasible.
   <h2>Index of all input sections</h2>
    <xsl:call-template name="section_index"></xsl:call-template>
    Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
   <hr/>
   <h2>Detailed description of all sections and keywords</h2>
    <xsl:call-template name="describe_sections"></xsl:call-template>
   <hr/>
   Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
  </body>
 </html>
</xsl:template>

<xsl:template name="section_index">
 <xsl:for-each select="SECTION">
  <xsl:sort select="NAME"/>
  <ul class="disc">
   <li><a href="#sec_des_{generate-id(NAME)}" id="sec_ind_{generate-id(NAME)}">&amp;<xsl:value-of select="NAME"/></a>
    <!--xsl:call-template name="keyword_index"></xsl:call-template-->
    <xsl:call-template name="section_index"></xsl:call-template>
   </li>
  </ul>
 </xsl:for-each>
</xsl:template>

<xsl:template name="keyword_index">
 <xsl:for-each select="KEYWORD">
  <xsl:sort select="NAME[@type='default']"/>
  <ul class="circle">
   <li><a href="#key_des_{generate-id(NAME[@type='default'])}" id="key_ind_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a></li>
  </ul>
 </xsl:for-each>
</xsl:template>

<xsl:template name="describe_sections">
 <xsl:for-each select="SECTION">
  <xsl:sort select="NAME"/>
  <p>
  <h3><a href="#sec_ind_{generate-id(NAME)}" id="sec_des_{generate-id(NAME)}">Section &amp;<xsl:value-of select="NAME"/></a></h3>
  <ul class="none">
   <li>
    <em><xsl:value-of select="DESCRIPTION"/></em>
    <xsl:if test="count(KEYWORD) > 0">
     <h4>Keywords:</h4>
     <xsl:call-template name="describe_keywords"></xsl:call-template>
    </xsl:if>
   </li>
  </ul>
  </p>
  <xsl:call-template name="describe_sections"></xsl:call-template>
 </xsl:for-each>
</xsl:template>

<xsl:template name="describe_keywords">
 <ul class="none">
  <li>
  <xsl:for-each select="KEYWORD">
   <xsl:sort select="NAME"/>
   <dl>
    <dt>
     <!--a href="#key_ind_{generate-id(NAME[@type='default'])}" id="key_des_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a-->
     <a href="#sec_ind_{generate-id(../NAME)}"><xsl:value-of select="NAME[@type='default']"/></a>
     <xsl:if test="NAME[@type='alias']">
      (alias:
      <xsl:for-each select="NAME[@type='alias']">
       <xsl:sort select="NAME[@type='alias']"/>
       <xsl:choose>
        <xsl:when test="position() = last()">
         <xsl:value-of select="."/>)
        </xsl:when>
        <xsl:otherwise>
         <xsl:value-of select="."/>,
        </xsl:otherwise>
       </xsl:choose>
      </xsl:for-each>
     </xsl:if>
    </dt>
    <dd><p><em><xsl:value-of select="DESCRIPTION"/></em></p></dd>
    <dd>
     <p>
      Data type: <big class="uctt"><xsl:value-of select="DATA_TYPE/@kind"/></big>
      <xsl:if test="DATA_TYPE/ENUMERATION">
       <p>List of valid keys:</p>
       <ul class="none">
       <xsl:for-each select="DATA_TYPE/ENUMERATION/ITEM">
        <xsl:sort select="NAME"/>
        <dl>
         <dt><big class="uctt"><xsl:value-of select="NAME"/></big></dt>
         <dd><em><xsl:value-of select="DESCRIPTION"/></em></dd>
        </dl>
       </xsl:for-each>
       </ul>
      </xsl:if>
     </p>
    </dd>
    <dd><p>Usage example: <big class="uctt"><xsl:value-of select="USAGE"/></big></p></dd>
    <xsl:if test="string-length(DEFAULT_VALUE) > 0">
     <dd><p>Default value: <big class="uctt"><xsl:value-of select="DEFAULT_VALUE"/></big></p></dd>
    </xsl:if>
   </dl>
  </xsl:for-each>
  </li>
 </ul>
</xsl:template>

</xsl:stylesheet>
