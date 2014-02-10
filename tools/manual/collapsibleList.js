/*****************************************************************************************************************
                Automatic collapsible list generation script written by Mark Wilton-Jones - 2002
Version 3.0 - updated May 2004 to allow cookies to remember the collapse state - fully compatible with original script
******************************************************************************************************************

Please see http://www.howtocreate.co.uk/jslibs/ for details and a demo of this script
Please see http://www.howtocreate.co.uk/jslibs/termsOfUse.html for terms of use
_________________________________________________________________________

You can put as lists on the page as you like, each list may have a different format

To use:
_________________________________________________________________________

Inbetween the <head> tags, put:

        <script src="PATH TO SCRIPT/collapsibleList.js" type="text/javascript" language="javascript1.2"></script>

If you want to be able to remember the expand/collapse state of all the branches so that it remains
open if you open a new page or return to the page, you will also need to include my cookie script
(available from http://www.howtocreate.co.uk/jslibs/)

        <script src="PATH TO SCRIPT/cookie.js" type="text/javascript" language="javascript1.2"></script>
_________________________________________________________________________

For best results, all images should be the same height, slightly taller than
the text around them. I would recommend about 21px.

Now, wherever you want the list (not in a layer or positioned element due
to Netscape 4 style bug), put:

        <script type="text/javascript" language="javascript1.2"><!--

//The following is an example of a list which I have decided to call 'myList'

//if using classes (optional), define the class array here (see below):
var classList = [
        ['classForLevel_0_WithNoChildBranch','classForLevel_0_WithChildBranch'],
        ['classForLevel_1_WithNoChildBranch','classForLevel_1_WithChildBranch'],
        ... etc ...
        ['classForLevel_n_WithNoChildBranch','classForLevel_n_WithChildBranch']
];

//create a tree object and define all images to use
var myList = new collapsibleList(
        [
                17,      //width of tree and branch pictures
                21,      //height of tree and branch pictures
                't.gif', //branch junction image, shaped like this |-
                'l.gif', //'last branch' image, shaped like this `-
                'i.gif', //tree trunk image, shaped like this |
                'e.gif'  //blank spacer image, eg. a transparent gif
        ],
        [
                25,      //width of expand / collapse images
                21,      //height of expand / collapse images
                'f.gif', //image used in place of expand or collapse images if the browser cannot expand / collapse
                'b.gif', //basic horizontal line image like this - for branches that do not expand
                'p.gif', //image used for expand link, looks like this [+]
                'm.gif'  //image used for collapse link, looks like this [-]
        ],
        false,     //optional - set to true to automatically collapse branches when sibling branches are
                   //expanded or parent branches are collapsed
        classList  //optional - array list of classes for branch levels (NOT applied in Netscape 4 because of browser bugs)
);

//then create branches
//format is:
//  myList.SUB_REF = new sub('String: HTML content'[,bool:expandByDefault]);
//expandByDefault is only used if list is not set to automatically collapse
myList.sub[0] = new sub('List item 1',true);
myList.sub[0].sub[0] = new sub('List <a href="canContainHTML.html">item</a> 1.1');
myList.sub[0].sub[0].sub[0] = new sub('List item 1.1.1');
myList.sub[0].sub[1] = new sub('List item 1.2',true);
myList.sub[0].sub[1].sub[0] = new sub('List item 1.2.1');
myList.sub[0].sub[1].sub[1] = new sub('List item 1.2.2',true);
myList.sub[0].sub[1].sub[1].sub[0] = new sub('List item 1.2.2.1');
myList.sub[0].sub[2] = new sub('List item 1.3');
myList.sub[0].sub[2].sub[0] = new sub('List item 1.3.1');
myList.sub[0].sub[3] = new sub('List item 1.4');
myList.sub[0].sub[4] = new sub('List item 1.5');
myList.sub[1] = new sub('List item 2');
myList.sub[2] = new sub('List item 3');
myList.sub[2].sub[0] = new sub('List item 3.1');
myList.sub[2].sub[0].sub[0] = new sub('List item 3.1.1');

//then tell the browser to create the collapsible list
createList(myList,'cookieNameToSave');

//'cookieNameToSave' is optional and will be used to save/recover the expand/collapse state of all branches
//State can only be saved for collapsible lists that are not set to automatically collapse branches when sibling
//branches are expanded or parent branches are collapsed

        //--></script>

_________________________________________________________________________

To trigger saving of expand/collapse state of all branches in all lists where you have specified a cookie name, use:

        saveCollapseState(lifeTime)

lifeTime is optional (defaults to browser session) and specifies the amount of time (in seconds)
to remember the expand/collapse state of all branches

Usually, this would be used onunload, eg. <body onunload="saveCollapseState();">

*******************************************************************************************************
                               And here's the actual code
******************************************************************************************************/

window.collapsibleListRef = []; window.imageCache = []; window.autoCloseRef = [];
//iCab can change display, but does not say so, and does notunderstand '' so I use a browser detect
var is_Icab = window.ScriptEngine && ScriptEngine().indexOf( 'InScript' ) + 1;

function collapsibleList() {
        //hold tree information about the list and the images used
        this.treeImg = arguments[0];
        this.expdImg = arguments[1];
        this.isAutoClose = arguments[2];
        this.subClass = arguments[3];
        this.sub = [];
}

function sub() {
        //hold information about the list item
        this.text = arguments[0];
        this.expanded = arguments[1];
        this.sub = [];
}

//done because WebTV does not understand inline object syntax
function ColLexPdOb(plI,miI,oSv,oAc) { this.plI = plI; this.miI = miI; this.subs = []; this.saveName = oSv; this.autoCol = oAc; }
function ColLexPdIT(name,expanded) { this.name = name; this.expanded = expanded; }

function createList(listObject,currentExt,imgExt,oLev,oBase) {
        if( !currentExt || !currentExt.match(/^treeList_\d+(_\d+)+$/) ) {
                //create the list, start by defining aditional stuff that is required, counting the lists etc.
                window.collapsibleListRef[window.collapsibleListRef.length] = new ColLexPdOb( listObject.expdImg[4], listObject.expdImg[5], currentExt, listObject.isAutoClose );
                currentExt = 'treeList_'+( window.collapsibleListRef.length - 1 ); imgExt = ''; oLev = 0; oBase = listObject;
                document.write( '<p style="line-height:'+(document.layers?listObject.treeImg[1]+'px':'100%')+';">' );
                oBase.menuUniqueNum = window.collapsibleListRef.length - 1; window.autoCloseRef[oBase.menuUniqueNum] = [];
        }
        document.write( '<span style="display:;" id="'+currentExt+'">' );
        for( var x = 0; x < listObject.sub.length; x++ ) {
                //for every child object, create their line
                document.write(
                        //the line break on all except the very first line
                        ( ( oLev || x ) ? '\n<br>' : '' ) +
                        //what may or may not remain of the trunk
                        ( ( oBase.subClass && oBase.subClass[oLev] && oBase.subClass[oLev][listObject.sub[x].sub.length?1:0] && !( document.layers && navigator.mimeTypes['*'] ) ) ? ( '<span class="'+oBase.subClass[oLev][listObject.sub[x].sub.length?1:0]+'">' ) : '' ) + imgExt +
                        //this branch
                        ( oLev ? '<img src="' + ( ( x < listObject.sub.length - 1 ) ? oBase.treeImg[2] : oBase.treeImg[3] ) + '" align="absmiddle" width="' + oBase.treeImg[0] + '" height="' + oBase.treeImg[1] + '" alt="' + ( ( x < listObject.sub.length - 1 ) ? '|-' : '`-' ) + '" border="0">' : '' ) +
                        //the expand / collapse link or not as the case may be
                        ( listObject.sub[x].sub.length ? '<a href="#" onclick="expandList(\'' + currentExt + '_' + x + '\',\'' + currentExt + '_' + x + '_img\',\'' + oBase.expdImg[4] + '\',\'' + oBase.expdImg[5] + '\',\'' + oBase.menuUniqueNum + '\','+oLev+',' + ( oBase.isAutoClose ? 'true' : 'false' ) + ');if(this.blur){this.blur();}return false;"><img src="' + oBase.expdImg[2] + '" align="absmiddle" name="' + currentExt + '_' + x + '_img" width="' + oBase.expdImg[0] + '" height="' + oBase.expdImg[1] + '" alt="[+/-]" border="0"></a> ' : '<img src="' + oBase.expdImg[3] + '" align="absmiddle" width="' + oBase.expdImg[0] + '" height="' + oBase.expdImg[1] + '" alt="-----" border="0"> ' ) +
                        //the text of the branch
                        listObject.sub[x].text + ( ( oBase.subClass && oBase.subClass[oLev] && oBase.subClass[oLev][listObject.sub[x].sub.length?1:0] && !( document.layers && navigator.mimeTypes['*'] ) ) ? '</span>' : '' )
                );
                if( listObject.sub[x].sub.length ) {
                        //add the span id to a list so we can easily find it and collapse it later
                        window.collapsibleListRef[window.collapsibleListRef.length - 1].subs[window.collapsibleListRef[window.collapsibleListRef.length - 1].subs.length] = new ColLexPdIT( currentExt + '_' + x, listObject.sub[x].expanded && !oBase.isAutoClose );
                        //create children
                        createList( listObject.sub[x], currentExt + '_' + x, oLev ? ( imgExt + ( ( x < listObject.sub.length - 1 ) ? '<img src="' + oBase.treeImg[4] + '" align="absmiddle" width="' + oBase.treeImg[0] + '" height="' + oBase.treeImg[1] + '" alt="|&nbsp;">' : '<img src="' + oBase.treeImg[5] + '" align="absmiddle" width="' + oBase.treeImg[0] + '" height="' + oBase.treeImg[1] + '" alt="&nbsp;&nbsp;">' ) ) : '', oLev + 1, oBase );
                }
        }
        document.write( '</span>' );
        if( !oLev ) {
                //end the list and prepare to collapse as soon as the browser lays it out
                document.write( '</p>\n' );
                if( document.all || document.getElementById ) { window.setTimeout('prepareForUse('+(window.collapsibleListRef.length-1)+')',100); }
        }
}

function expandList(spanName,imgName,plsImg,minImg,oMenu,oLevel,oAutoClose) {
        var theSpan = document.getElementById ? document.getElementById(spanName) : document.all ? document.all[spanName] : (new Object());
        if( !theSpan ) { return; } if( theSpan.style ) { theSpan = theSpan.style; }
        if( typeof( theSpan.display ) == 'undefined' && !is_Icab ) {
                //if we could not access the span element or if its display style cannot be changed, say so
                window.alert( 'Sorry, your browser cannot collapse or expand these lists\nso the lists will remain constantly expanded.' );
        } else {
                //if required, collapse back down again
                for( var y = window.autoCloseRef[oMenu].length - 1; y >= oLevel; y-- ) {
                        if( window.autoCloseRef[oMenu][y] && oAutoClose && window.autoCloseRef[oMenu][y][0] != spanName ) {
                                var theSpan2 = document.getElementById ? document.getElementById(window.autoCloseRef[oMenu][y][0]) : document.all ? document.all[window.autoCloseRef[oMenu][y][0]] : (new Object());
                                if( theSpan2.style ) { theSpan2 = theSpan2.style; } theSpan2.display = 'none';
                                document.images[window.autoCloseRef[oMenu][y][1]].src = window.autoCloseRef[oMenu][y][2];
                                window.autoCloseRef[oMenu][y] = null;
                        }
                }
                if( theSpan.display == 'none' ) { window.autoCloseRef[oMenu][oLevel] = [spanName,imgName,plsImg]; }
                //expand / collapse the list and change the image
                theSpan.display = ( theSpan.display == 'none' ) ? ( is_Icab ? 'inline' : '' ) : 'none';
                document.images[imgName].src = theSpan.display ? plsImg : minImg;
        }
}

function prepareForUse(listNum) {
        if( !window.collapsibleListRef[listNum].subs.length ) { return; } //no branches to collapse
        var lastPart = window.collapsibleListRef[listNum].subs[window.collapsibleListRef[listNum].subs.length - 1].name;
        lastPart = document.getElementById ? document.getElementById(lastPart) : document.all[lastPart];
        //if the page has not loaded enough, try again later
        if( !lastPart ) { window.setTimeout('prepareForUse('+listNum+')',100); } else {
                //if their browser cannot change the display style, don't try
                if( lastPart.style ) { lastPart = lastPart.style; }
                if( typeof( lastPart.display ) == 'undefined' && !is_Icab ) { return; }
                //cache images for faster changes
                window.imageCache[listNum] = [(new Image()),(new Image())];
                window.imageCache[listNum][0].src = window.collapsibleListRef[listNum].plI;
                window.imageCache[listNum][1].src = window.collapsibleListRef[listNum].miI;
                var svNm; if( svNm = window.collapsibleListRef[listNum].saveName ) { svNm = retrieveCookie(svNm); }
                if( typeof( svNm ) == 'string' ) { svNm = svNm.split(','); }
                for( var x = 0; x < window.collapsibleListRef[listNum].subs.length; x++ ) {
                        var lastPart = window.collapsibleListRef[listNum].subs[x];
                        if( svNm ) { var thisArray = mwjInArray(lastPart.name,svNm); }
                        if( ( svNm && thisArray ) || ( !svNm && lastPart.expanded ) ) {
                                //if they want it expanded by default, just change the image
                                document.images[lastPart.name + '_img'].src = window.imageCache[listNum][1].src;
                        } else {
                                //collapse the branch and change the image
                                document.images[lastPart.name + '_img'].src = window.imageCache[listNum][0].src;
                                lastPart = document.getElementById ? document.getElementById(lastPart.name) : document.all[lastPart.name];
                                if( lastPart.style ) { lastPart = lastPart.style; }
                                lastPart.display = 'none';
                        }
                }
        }
}

function mwjInArray(oNeed,oHay) { for( var i = 0; i < oHay.length; i++ ) { if( oNeed == oHay[i] ) { return true; } } return false; }

function saveCollapseState(oLife) {
        var allSpans = document.getElementsByTagName ? document.getElementsByTagName('span') : ( document.all ? document.all.tags('SPAN') : [] );
        if( !allSpans.length ) { return; } //nothing to save
        //try to minimise the number of comparisons, to save processor time
        for( var x = 0, saveNums = ''; window.collapsibleListRef[x]; x++ ) {
                if( !window.collapsibleListRef[x].autoCol && window.collapsibleListRef[x].saveName ) { saveNums += ( saveNums ? '|' : '' ) + x; } }
        if( !saveNums ) { return; } //nothing to save
        for( var x = 0, brnch = [], frNm, rgxp = new RegExp('^treeList_('+saveNums+')_.*$',''); x < allSpans.length; x++ ) {
                if( allSpans[x].id && allSpans[x].id.match(rgxp) && allSpans[x].style.display != 'none' ) {
                        frNm = allSpans[x].id.replace(rgxp,'$1');
                        brnch[frNm] = ( brnch[frNm] ? ( brnch[frNm] + ',' ) : '' ) + allSpans[x].id;
                }
        }
        for( var x = 0, saveNums = ''; window.collapsibleListRef[x]; x++ ) {
                if( window.collapsibleListRef[x].saveName ) { setCookie(window.collapsibleListRef[x].saveName,brnch[x],oLife,'/'); } }
}
