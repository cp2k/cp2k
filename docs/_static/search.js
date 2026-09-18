if (window.location.protocol !== "file:") {
    const root = new URL("../", document.currentScript.src);
    const bundle = new URL("pagefind/", root);
    const config = document.querySelector("pagefind-config");

    config?.setAttribute("base-url", root.pathname);
    config?.setAttribute("bundle-path", bundle.pathname);

    customElements.whenDefined("pagefind-modal-trigger").then(() => {
        document.getElementById("cp2k-pagefind-search").hidden = false;
        document.getElementById("cp2k-classic-search").hidden = true;
    });
}
