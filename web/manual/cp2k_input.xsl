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
   <script src="collapsibleList.js" type="text/javascript" language="javascript1.2">
    this line is here just to avoid that Xalan merges the previous and the next line
   </script>
  </head>
  <body>
   <h1 align="center">CP2K input reference</h1>
   <h2>Version information</h2>
   <p>
    This HTML manual was generated automatically from a CP2K executable
    compiled on <xsl:value-of select="COMPILE_DATE"/> using the
    --xml command line option. Thus the manual describes exactly this
    version of the code. The latest CVS log file entry found was
   <xsl:value-of select="COMPILE_LASTCVS"/>.</p>

   <h2>Journal papers</h2>
   <p><a href="references.html">List of references</a> cited in the CP2K input manual.</p>

   <h2>Internal Input Preprocessor</h2>
   <p>Before the input is parsed, the input is run through a simple internal preprocessor.
   The preprocessor recognizes the following directives independent of capitalization:</p>
   <dl>
     <dt><b>@INCLUDE 'filename.inc'</b></dt>
     <DD>The file referenced by <i>filename.inc</i> is included into the input file and parsed.
     Recursive inclusions are not allowed and the files have to exist in the current working 
     directory. There can be only one @INCLUDE statement per line. Single or double quotes 
     have to be used if the filename contains blanks.</DD> 
     <DT><b>@SET VAR value</b></DT>
     <DD>Assigns the text <i>value</i> to the preprocessing variable <i>VAR</i>. <i>value</i>
     is the text following <i>VAR</i> with the outer whitespace removed. The variable can be
     recalled with a <i>${VAR}</i> statement. There can be only one @SET statement per line.</DD>
     <DT><b>${VAR}</b></DT>
     <DD>Expand the variable <i>VAR</i>. The text <i>${VAR}</i> is replaced with the value assigned 
     to <i>VAR</i> in the last @SET directive. There can be multiple variable statements
     per line. The expansion process is repeated until no more variables are found.</DD>
     <DT><B>@IF / @ENDIF</B></DT>
     <DD>Conditional block. The text from the @IF line up to the next line with a valid "
     @ENDIF is skipped, if the expression following @IF resolves to <I>false</I>. 
     Available expressions are lexical comparisons for equality '==' or inequality '/='.
     If none of the two operators are found, a '0' or whitespace resolves to <I>false</I>
     while any text resolves to <I>true</I>. @IF/@ENDIF blocks cannot be nested and
     cannot span across files. There can be only one test (== or /=) per @IF statement.</DD>
   </dl>

   <h2>Input structure</h2>
    All sections and keywords that can be part of a CP2K input file are shown
    with their allowed nestings. Click on <img src="p.gif" alt="expand"/> to
    expand and on <img src="m.gif" alt="collapse"/> to collapse the current
    index level. A detailed description can be obtained by clicking on the
    links. The links in the detailed descriptions switch back to the
    corresponding index entries if they are currently expanded. In this way
    a toggling between the index and the detailed description is feasible.
   <h2>Index of all input sections and keywords</h2>
    <h4>Last update: <xsl:value-of select="COMPILE_DATE"/></h4>
    <script type="text/javascript" language="javascript1.2">
     var Path = new collapsibleList([17,21,'t.gif','l.gif','i.gif','e.gif'],[25,21,'f.gif','b.gif','p.gif','m.gif'],false);
     <xsl:call-template name="section_index">
      <xsl:with-param name="path" select="'Path'"/>
      <xsl:with-param name="expand_list" select="'false'"/>
      <xsl:with-param name="start" select="count(KEYWORD)"/>
     </xsl:call-template>
     createList(Path);
    </script>
    Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
   <hr/>
   <h2>Detailed description of all sections and keywords</h2>
    <h4>Last update: <xsl:value-of select="COMPILE_DATE"/></h4>
    <xsl:call-template name="describe_sections">
     <xsl:with-param name="path" select="''"/>
    </xsl:call-template>
   <hr/>
   Back to the <a href="http://cp2k.berlios.de">CP2K home page</a>
  </body>
 </html>
</xsl:template>

<xsl:template name="section_index">
 <xsl:param name="path"/>
 <xsl:param name="expand_list"/>
 <xsl:param name="start"/>
 <xsl:for-each select="SECTION">
  <xsl:sort select="NAME"/>
  <xsl:variable name="local_path" select="concat($path,'.sub[',string($start+position()-1),']')"/>
  <xsl:value-of select="$local_path"/> = new sub('<a href="#sec_des_{generate-id(NAME)}" id="sec_ind_{generate-id(NAME)}">&amp;<xsl:value-of select="NAME"/></a>',<xsl:value-of select="$expand_list"/>);
  <xsl:call-template name="keyword_index">
   <xsl:with-param name="path" select="$local_path"/>
  </xsl:call-template>
  <xsl:call-template name="section_index">
   <xsl:with-param name="path" select="$local_path"/>
   <xsl:with-param name="expand_list" select="$expand_list"/>
   <xsl:with-param name="start" select="count(SECTION_PARAMETERS)+count(DEFAULT_KEYWORD)+count(KEYWORD)"/>
  </xsl:call-template>
 </xsl:for-each>
</xsl:template>

<xsl:template name="keyword_index">
 <xsl:param name="path"/>
 <xsl:variable name="nsecpar" select="count(SECTION_PARAMETERS)"/>
 <xsl:variable name="ndefpar" select="count(DEFAULT_KEYWORD)"/>
 <xsl:for-each select="SECTION_PARAMETERS">
  <xsl:sort select="NAME[@type='default']"/>
  <xsl:variable name="local_path" select="concat($path,'.sub[0]')"/>
  <xsl:value-of select="$local_path"/> = new sub('<a href="#key_des_{generate-id(NAME[@type='default'])}" id="key_ind_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a>');
 </xsl:for-each>
 <xsl:for-each select="DEFAULT_KEYWORD">
  <xsl:sort select="NAME[@type='default']"/>
  <xsl:variable name="local_path" select="concat($path,'.sub[',string($nsecpar),']')"/>
  <xsl:value-of select="$local_path"/> = new sub('<a href="#key_des_{generate-id(NAME[@type='default'])}" id="key_ind_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a>');
 </xsl:for-each>
 <xsl:for-each select="KEYWORD">
  <xsl:sort select="NAME[@type='default']"/>
  <xsl:variable name="local_path" select="concat($path,'.sub[',string($nsecpar+$ndefpar+position()-1),']')"/>
  <xsl:value-of select="$local_path"/> = new sub('<a href="#key_des_{generate-id(NAME[@type='default'])}" id="key_ind_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a>');
 </xsl:for-each>
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
    <big class="uctt">
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
  <!--ul class="none">
   <li>
    This section is required:
    <big class="uctt">
     <xsl:value-of select="@required"/>
    </big>
   </li>
  </ul-->
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
   <ul class="none">
    <li>
     <h3>Keywords:</h3>
     <xsl:call-template name="describe_keywords"/>
    </li>
   </ul>
  </xsl:if>
  <xsl:call-template name="describe_sections">
   <xsl:with-param name="path" select="$local_path"/>
  </xsl:call-template>
 </xsl:for-each>
</xsl:template>

<xsl:template name="describe_keywords">
 <xsl:param name="element" select="KEYWORD"/>
 <ul class="none">
  <li>
   <xsl:for-each select="$element">
    <xsl:sort select="NAME"/>
    <xsl:if test="not(starts-with(NAME,'__'))">
     <dl>
      <dt>
       <a href="#key_ind_{generate-id(NAME[@type='default'])}" id="key_des_{generate-id(NAME[@type='default'])}"><xsl:value-of select="NAME[@type='default']"/></a>
       <!--a href="#sec_ind_{generate-id(../NAME)}"><xsl:value-of select="NAME[@type='default']"/></a-->
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
      <!--dd>
       This keyword is required:
       <big class="uctt">
        <xsl:value-of select="@required"/>
       </big>
      </dd-->
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
   </xsl:for-each>
  </li>
 </ul>
</xsl:template>

</xsl:stylesheet>
