<?xml version="1.0" encoding="ISO-8859-1"?>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<xsl:template match="/INPUT/SECTION">
 <html>
  <head>
   <title>CP2K input</title>
  </head>
  <body bgcolor="#e0e0e0">
   <h1 align="center">CP2K input description</h1>
   <hr/>
   <h2>List of sections and keywords</h2>
   <xsl:call-template name="section_index"></xsl:call-template>
   <hr/>
   <h2>Detailed description of all sections and keywords</h2>
   <xsl:call-template name="describe_sections"></xsl:call-template>
   <hr/>
  </body>
 </html>
</xsl:template>

<xsl:template name="section_index">
 <xsl:choose>
  <xsl:when test="count(SECTION) = 0">
   <br/>
  </xsl:when>
  <xsl:otherwise>
   <xsl:for-each select="SECTION">
    <xsl:sort select="NAME"/>
    <ul type="disc">
     <li><a href="#sec_des_{generate-id(NAME)}" id="sec_ind_{generate-id(NAME)}">&amp;<xsl:value-of select="NAME"/></a></li>
     <xsl:call-template name="keyword_index"></xsl:call-template>
     <xsl:call-template name="section_index"></xsl:call-template>
    </ul>
   </xsl:for-each>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<xsl:template name="keyword_index">
 <xsl:for-each select="KEYWORD">
  <xsl:sort select="NAME[@type='default']"/>
  <ul type="circle">
   <li><a href="#key_des_{generate-id(NAME[@type='default'])}" id="key_ind_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a></li>
  </ul>
 </xsl:for-each>
 <br/>
</xsl:template>

<xsl:template name="describe_sections">
 <xsl:choose>
  <xsl:when test="count(SECTION) = 0">
  </xsl:when>
  <xsl:otherwise>
   <xsl:for-each select="SECTION">
    <xsl:sort select="NAME"/>
    <br/>
    <a href="#sec_ind_{generate-id(NAME)}" id="sec_des_{generate-id(NAME)}"><h3>Section &amp;<xsl:value-of select="NAME"/></h3></a>
    <ul>
     <xsl:value-of select="DESCRIPTION"/>
     <xsl:call-template name="describe_keywords"></xsl:call-template>
    </ul>
    <xsl:call-template name="describe_sections"></xsl:call-template>
   </xsl:for-each>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

<xsl:template name="describe_keywords">
 <xsl:choose>
  <xsl:when test="count(KEYWORD) = 0">
  </xsl:when>
  <xsl:otherwise>
   <h4>Keywords:</h4>
   <ul>
   <xsl:for-each select="KEYWORD">
    <xsl:sort select="NAME"/>
    <dl>
     <xsl:if test="position() != 1"><br/></xsl:if>
     <dt>
      <a href="#key_ind_{generate-id(NAME[@type='default'])}" id="key_des_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a>
      <xsl:if test="NAME[@type='alias']">
       (alias: <xsl:value-of select="NAME[@type='alias']"/>)
      </xsl:if>
     </dt> 
     <dd><p><em><xsl:value-of select="DESCRIPTION"/></em></p></dd>
     <dd>
      <p>Data type:
      <xsl:choose>
       <xsl:when test="DATA_TYPE/BOOLEAN">
        <big><tt>BOOLEAN</tt></big>
       </xsl:when>
       <xsl:when test="DATA_TYPE/INTEGER">
        <big><tt>INTEGER</tt></big>
       </xsl:when>
       <xsl:when test="DATA_TYPE/FLOAT">
        <big><tt>FLOAT</tt></big>
       </xsl:when>
       <xsl:when test="DATA_TYPE/STRING">
        <big><tt>STRING</tt></big>
       </xsl:when>
       <xsl:when test="DATA_TYPE/ENUMERATION">
        <big><tt>ENUMERATION</tt></big>
        <p>List of valid keys:<br/>
        <ul>
         <xsl:for-each select="DATA_TYPE/ENUMERATION/ITEM">
          <xsl:sort select="NAME"/>
          <dl>
           <dt><big><tt><xsl:value-of select="NAME"/></tt></big></dt>
           <dd><xsl:value-of select="DESCRIPTION"/></dd>
          </dl>
         </xsl:for-each>
        </ul>
        </p>
       </xsl:when>
       <xsl:otherwise>
        <big><tt>undefined</tt></big>
       </xsl:otherwise>
      </xsl:choose>
      </p>
     </dd>
     <dd><p>Usage example: <big><tt><xsl:value-of select="USAGE"/></tt></big></p></dd>
     <xsl:if test="string-length(DEFAULT_VALUE) > 0">
      <dd><p>Default value: <big><tt><xsl:value-of select="DEFAULT_VALUE"/></tt></big></p></dd>
     </xsl:if>
    </dl>
   </xsl:for-each>
   </ul>
  </xsl:otherwise>
 </xsl:choose>
</xsl:template>

</xsl:stylesheet>
