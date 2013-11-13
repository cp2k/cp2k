/****** toggle_folding/main
 * PURPOSE
 *   Allow collapsible section tree in CP2K manual
 * USAGE
 *   load this script in the header of the CP2K manual font page or
 *   the template file.
 *
 *   If javascript is enabled, then the user should see [-] buttons in
 *   front of all section titles. And:
 *
 *   - Single clicking on a [-]/[+] button folds/unfolds the
 *     subsection tree
 *   - Double clicking on a [-]/[+] button folds/unfolds the entire
 *     subsection tree recursively.
 *
 *  The buttons are defined as HTML <a> elements, and have class name
 *  "button", and one can change its font and appearance by setting
 *  the corresponding style for .button in the CSS.
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
// add the folding buttons to the manual tree
document.addEventListener('DOMContentLoaded', parse_manual, false);
/*****/

/****** toggle_folding/parse_manual
 * PURPOSE
 *   Parses the CP2K manual front page and addes the tree folding
 *   buttons to section titles
 * USAGE
 *   parse_manual();
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function parse_manual() {
    // find the first element with id "CP2K_INPUT.html"
    var root = document.getElementById("CP2K_INPUT.html");
    // the parent element of root is the section title element
    var root_section = root.parentNode;
    // parse the root_section recursively
    add_button_recursively(root_section);
    return;
}
/*****/

/****f* toggle_folding/add_button_recursively
 * PURPOSE
 *   Recursively adds the tree-folding buttons to a section tree
 *   Inserts a blank button if a section does not subsections
 * USAGE
 *   add_button_recursively(section)
 * INPUTS
 *   section: the HTML <li> element in the CP2K input list containing
 *            the section title
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function add_button_recursively(section) {
    var button;
    var title_element = first_element_child(section);
    var subsections = next_element_sibling(section);
    var subsection = subsections ? first_element_child(subsections) : null;
    // only add button if the subsections are not empty
    if (title_element) {
        // add button if the subsections are not empty
        if (subsection) {
            button = make_button("[\u2212]", "toggle");
            section.insertBefore(button, title_element);
        }
        // otherwise add a blank
        else {
            button = make_button("\u00a0\u00a0\u00a0", "blank");
            section.insertBefore(button, title_element);
        }
    }
    // do this recursively for all subsections
    while (subsection) {
        // the element is a section title and only needs parsing if it
        // contains an link element with a non-null id
        if (first_element_child(subsection).id) {
            add_button_recursively(subsection);
            // skip the subsections of the subsection
            subsection = next_element_sibling(subsection);
        }
        subsection = next_element_sibling(subsection);
    }
    return;
}
/*****/

/****f* toggle_folding/make_button
 * PURPOSE
 *   Makes a button element, in the form of
 *   <a onclick="toggle(event)" data-dblclick="nil" data-folding="open">text</a>
 *
 *   If the button is of the type "blank", then it will not be
 *   clickable, and will have the data-folding set to blank. The blank
 *   buttons are used as place-holders to keep of section titles
 *   properly aligned.
 * USAGE
 *   make_button("[-]")
 * INPUTS
 *   text: the text "image" of the button
 *   type: "toggle" | "blank"
 * RETURN VALUE
 *   a new <a> element representing a button
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function make_button(text, type) {
    var button;
    var button_text;
    button = document.createElement("a");
    button_text = document.createTextNode(text.concat(" "));
    button.appendChild(button_text);
    button.className = "button";
    switch (type) {
    case "toggle":
        // add onclick attribute
        button.setAttribute("onclick", "toggle(event)");
        // for recording double click
        button.setAttribute("data-dblclick", "nil");
        // for recording the open and closed status
        button.setAttribute("data-folding", "open");
        break;
    case "blank":
    default :
        button.setAttribute("data-folding", "blank");
    }
    return button;
}
/*****/

/****f* toggle_folding/toggle_button
 * PURPOSE
 *   Toggle the button element to open or close status
 *   Do nothing if the button is blank
 * USAGE
 *   toggle_button(button, "open")
 * INPUTS
 *   button: the <a> element representing button to toggle
 *   openOrClose:  "open" | "close"
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function toggle_button(button, openOrClose) {
    if (button.getAttribute("data-folding") === "blank")
        return;
    if (openOrClose === "open") {
        button.setAttribute("data-folding", "open");
        button.innerHTML = button.innerHTML.replace('[+]', '[\u2212]');
    }
    else if (openOrClose === "close") {
        button.setAttribute("data-folding", "close");
        button.innerHTML = button.innerHTML.replace('[\u2212]', '[+]');
    }
    return;
}
/*****/

/****f* toggle_folding/toggle
 * PURPOSE
 *   Event interface function for putting functionality for folding and
 *   unfolding the section trees in CP2K manual.
 *
 *   - single click on the [-] or [+] in front of section headings
 *     should fold/unfold the subtree corresponding to the section.
 *
 *   - double click on the [-] or [+] in front of section headings
 *     should fold/unfold the section subtree recursively.
 *
 *   Note that one cannot use ondbclick attribute in HTML togeter with
 *   onclick, as onclick will always fire first upon the first mouse
 *   click. This function provides a work-around solution to this
 *   problem.
 *
 *   The local variable delay controls the interval between the
 *   double-clicks.
 * USAGE
 *   Set onclick attribute to toggle(event)
 * INPUTS
 *   Mouse clicking event
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function toggle(event) {
    var delay = 300; // in milliseconds
    var el = event.target;
    var section = el.parentNode;
    if (el.getAttribute("data-dblclick") === "nil") {
        el.setAttribute("data-dblclick", "t");
        setTimeout(
            function() {
                if (el.getAttribute("data-dblclick") === "t") {
                    // single click registered
                    toggle_folding(section, false);
                }
                el.setAttribute("data-dblclick", "nil");
            },
            delay);
    }
    else {
        el.setAttribute("data-dblclick", "nil");
        // double click registered
        toggle_folding(section, true);
    }
    return;
}
/*****/

/****f* toggle_folding/toggle_folding
 * PURPOSE
 *   Toggling folding and unfolding of subsections, recursively if
 *   asked to.
 * USAGE
 *   toggle_folding(section, recursive)
 * INPUTS
 *   section: the HTML <li> element containing the section title,
 *            corresponds to the section to be folded/unfolded
 *   recursive:  true | false,  whether to do it recursively or not
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function toggle_folding(section, recursive) {
    if (! section)
        return;
    // get the next non-text sibling
    var subsections = next_element_sibling(section);
    if (! subsections)
        return;
    var button = first_element_child(section)
    if (button.className !== "button")
        return;
    if (button.getAttribute("data-folding") === "close")
        toggle_subtree(section, "open", recursive)
    else if (button.getAttribute("data-folding") === "open")
        toggle_subtree(section, "close", recursive)
    return;
}
/*****/

/****f* toggle_folding/toggle_subtree
 * PURPOSE
 *   The main driver function for doing the folding and unfolding of
 *   subsections
 *
 *   If openOrClose == "open", then open the section tree, and
 *   recursively open the section trees of the subsections if
 *   recursive == true. Similarly, openOrClose == "close" closes the
 *   section tree and recursively so if recursive == true.
 *
 *   This function assumes that: the non-text element where the
 *   folding function is triggered is always the title of the section.
 *   The next non-text element should be the list of subsections. All
 *   one needs to do, therefore, is to hide the subsections, i.e. the
 *   next non-text element.
 * USAGE
 *   toggle_folding(section, "open", false)
 * INPUTS
 *   section: the HTML <li> element containing the section title,
 *            corresponds to the section to be folded/unfolded
 *   openOrClose: "open" | "close"
 *   recursive:   true | false,  whether to do it recursively or not
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function toggle_subtree(section, openOrClose, recursive) {
    if (! section)
        return;
    var subsections = next_element_sibling(section);
    if (! subsections)
        return;
    var button = first_element_child(section)
    if (button.className !== "button")
        return;
    if (openOrClose === "open") {
        subsections.style.display = '';
        toggle_button(button, "open");
    }
    else if (openOrClose === "close") {
        subsections.style.display = "none";
        toggle_button(button, "close");
    }
    else {
        return;
    }
    // (un)fold recursively
    if (recursive) {
        var subsection = first_element_child(subsections);
        while (subsection) {
            // dump(subsection)
            if (first_element_child(subsection) &&
                (first_element_child(subsection).className === "button")) {
                toggle_subtree(subsection, openOrClose, true);
                subsection = next_element_sibling(subsection);
            }
            subsection = next_element_sibling(subsection);
        }
    }
    return;
}
/*****/

/****f* toggle_folding/first_element_child
 * PURPOSE
 *   find the first non-text (element) child of a HTML element
 * USAGE
 *   child = first_element_child(element)
 * INPUTS
 *   node: the HTML element in question
 * RETURN VALUE
 *   the first non-text child of node
 *   returns null if node has no children
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function first_element_child(node) {
    var el = node.firstChild;
    while (el && el.nodeType !== 1)
        el = el.nextSibling;
    return el;
}
/*****/

/****f* toggle_folding/next_element_sibling
 * PURPOSE
 *   find the next non-text (element) sibling of a HTML element
 * USAGE
 *   sibling = next_element_sibling(element)
 * INPUTS
 *   node: the HTML element in question
 * RETURN VALUE
 *   the next non-text (element) sibling of node
 *   returns null if node is the last element of its parent
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function next_element_sibling(node) {
    var el = node.nextSibling;
    while (el && el.nodeType !== 1)
        el = el.nextSibling;
    return el;
}
/*****/

/****f* toggle_folding/next_element_sibling
 * PURPOSE
 *   For debugging. Dumps the contents of a HTML element via alert
 * USAGE
 *   dump(element)
 * INPUTS
 *   element: the HTML element in question
 * AUTHOR
 *   L.Tong
 * CREATION DATE
 *   2013/11/06
 * MODIFICATION HISTORY
 * SOURCE
 */
function dump(element) {
  var a = ["Element dump:"];
  for (var k in element) {
    if (element.hasOwnProperty(k)) {
      a.push(k + ": " + element[k]);
    }
  }
  a.push("HTML: " + element.innerHTML);
  alert(a.join('\n'));
}
/*****/
