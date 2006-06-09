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
    body {background-color: #ffffff}
    big.tt {font-family: monospace}
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
    All sections and keywords that can be part of a CP2K input file are shown
    with their allowed nestings. A detailed description can be obtained by
    clicking on the section links. The links in the detailed descriptions switch
    back to the corresponding index entries. In this way a toggling between the
    index and the detailed description is feasible.
   <h2>Index of all input sections</h2>
    <h4>Last update: <xsl:value-of select="COMPILE_DATE"/></h4>
    <xsl:call-template name="section_index"/>
    Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
   <hr/>
   <h2>Detailed description of all sections and keywords</h2>
    <h4>Last update: <xsl:value-of select="COMPILE_DATE"/></h4>
    <xsl:call-template name="describe_sections"/>
  </body>
 </html>
</xsl:template>

<xsl:template name="section_index">
 <ul class="disc">
  <xsl:for-each select="SECTION">
   <xsl:sort select="NAME"/>
   <li>
    <a href="#sec_des_{generate-id(NAME)}" id="sec_ind_{generate-id(NAME)}"><xsl:value-of select="NAME"/></a>
   </li>
   <xsl:call-template name="section_index"/>
  </xsl:for-each>
 </ul>
</xsl:template>

<xsl:template name="subsection_index">
 <ul class="disc">
  <xsl:for-each select="SECTION">
   <xsl:sort select="NAME"/>
   <li>
    <a href="#sec_des_{generate-id(NAME)}" id="sec_ind_{generate-id(NAME)}"><xsl:value-of select="NAME"/></a>
   </li>
  </xsl:for-each>
 </ul>
</xsl:template>

<xsl:template name="describe_sections">
 <xsl:param name="path"/>
 <xsl:for-each select="SECTION">
  <xsl:sort select="NAME"/>
  <xsl:variable name="local_path" select="concat($path,'/',string(NAME))"/>
  <h3><a href="#sec_ind_{generate-id(NAME)}" id="sec_des_{generate-id(NAME)}">Section &amp;<xsl:value-of select="NAME"/></a></h3>
  <xsl:if test="string-length(DESCRIPTION) > 0">
   <ul class="none">
    <li>
     <em>
      <xsl:value-of select="DESCRIPTION"/>
     </em>
    </li>
   </ul>
  </xsl:if>
  <ul class="none">
   <li>
    Section path:
    <big class="tt">
     <xsl:value-of select="$local_path"/>
    </big>
   </li>
  </ul>
  <ul class="none">
   <li>
    This section can be repeated:
    <big class="uctt">
     <xsl:value-of select="@repeats"/>
    </big>
   </li>
  </ul>
  <ul class="none">
   <li>
    This section is required:
    <big class="uctt">
     <xsl:value-of select="@required"/>
    </big>
   </li>
  </ul>
  <ul class="none">
   <li>
    <h3>Subsections</h3>
    <xsl:choose>
     <xsl:when test="count(SECTION) > 0">
      <xsl:call-template name="subsection_index"/>
     </xsl:when>
     <xsl:otherwise>
      <ul class="none">
       <li>
        none
       </li>
      </ul>
     </xsl:otherwise>
    </xsl:choose>
   </li>
  </ul>
  <ul class="none">
   <li>
    <h3>Keywords:</h3>
    <xsl:choose>
     <xsl:when test="count(SECTION_PARAMETERS) > 0 or count(DEFAULT_KEYWORD) > 0 or count(KEYWORD) > 0">
      <xsl:if test="count(SECTION_PARAMETERS) > 0">
       <xsl:call-template name="describe_keywords">
        <xsl:with-param name="element" select="SECTION_PARAMETERS"/>
       </xsl:call-template>
      </xsl:if>
      <xsl:if test="count(DEFAULT_KEYWORD) > 0">
       <xsl:call-template name="describe_keywords">
        <xsl:with-param name="element" select="DEFAULT_KEYWORD"/>
       </xsl:call-template>
      </xsl:if>
      <xsl:if test="count(KEYWORD) > 0">
       <xsl:call-template name="describe_keywords"/>
      </xsl:if>
     </xsl:when>
     <xsl:otherwise>
      <ul class="none">
       <li>
        none
       </li>
      </ul>
     </xsl:otherwise>
    </xsl:choose>
   </li>
  </ul>
  Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
  <hr/>
  <xsl:call-template name="describe_sections">
   <xsl:with-param name="path" select="$local_path"/>
  </xsl:call-template>
 </xsl:for-each>
</xsl:template>

<xsl:template name="describe_keywords">
 <xsl:param name="element" select="KEYWORD"/>
 <ul class="none">
  <xsl:for-each select="$element">
   <xsl:sort select="NAME"/>
   <li>
    <xsl:if test="not(starts-with(NAME,'__'))">
     <dl>
      <dt>
       <u><xsl:value-of select="NAME[@type='default']"/></u>
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
      <dd>
       <p>
        <em>
         <xsl:value-of select="DESCRIPTION"/>
        </em>
       </p>
      </dd>
      <dd>
       This keyword can be repeated:
       <big class="uctt">
        <xsl:value-of select="@repeats"/>
       </big>
      </dd>
      <dd>
       <p>
       This keyword is required:
       <big class="uctt">
        <xsl:value-of select="@required"/>
       </big>
       </p>
      </dd>
      <dd>
       <p>
        Data type: <big class="uctt"><xsl:value-of select="DATA_TYPE/@kind"/></big>
       </p>
       <xsl:if test="DATA_TYPE/ENUMERATION">
        <p>
         List of valid keys:
        </p>
        <ul class="none">
         <li>
          <xsl:for-each select="DATA_TYPE/ENUMERATION/ITEM">
           <xsl:sort select="NAME"/>
           <dl>
            <dt>
             <big class="uctt">
              <xsl:value-of select="NAME"/>
             </big>
            </dt>
            <xsl:if test="string-length(DESCRIPTION) > 0">
             <dd>
              <em>
               <xsl:value-of select="DESCRIPTION"/>
              </em>
             </dd>
            </xsl:if>
           </dl>
          </xsl:for-each>
         </li>
        </ul>
       </xsl:if>
      </dd>
      <xsl:if test="string-length(USAGE) > 0">
       <dd>
        <p>
         Usage example:
         <big class="uctt">
          <xsl:value-of select="USAGE"/>
         </big>
        </p>
       </dd>
      </xsl:if>
      <xsl:if test="string-length(DEFAULT_VALUE) > 0">
       <dd>
        <p>
         Default value:
         <big class="uctt">
          <xsl:value-of select="DEFAULT_VALUE"/>
         </big>
        </p>
       </dd>
      </xsl:if>
      <xsl:if test="string-length(DEFAULT_UNIT) > 0">
       <dd>
        <p>
         Default unit:
         <big class="uctt">
          <xsl:value-of select="DEFAULT_UNIT"/>
         </big>
        </p>
       </dd>
      </xsl:if>
     </dl>
    </xsl:if>
   </li>
  </xsl:for-each>
 </ul>
</xsl:template>

</xsl:stylesheet>
