
window.ModelExporter = {
    name: "model",
    assertModelNonempty: function() {
        let isEmpty = Module.ccall('isModelEmpty', 'int', [], []);
        if (isEmpty) {
            alert("Model is empty.");
            return true;
        }
        return false;
    },
    downloadFile: function(ptr, filename, type) {
        let size = Module.ccall('getFileSize', 'int', [], []);
        let resultArray = new Uint8Array(Module.HEAPU8.buffer, ptr, size);
        let resultBuffer = resultArray.slice().buffer;
        console.log(size + " bytes");
        let link = document.createElement('a');
        link.href = URL.createObjectURL(new Blob([resultBuffer], { type: type }));
        link.download = filename;
        link.click();
        URL.revokeObjectURL(link.href);
    },
    downloadSTL: function() {
        if (ModelExporter.assertModelNonempty()) return;
        let ptr = Module.ccall('generateSTL', 'int', [], []);
        ModelExporter.downloadFile(ptr, ModelExporter.name+'.stl', 'model/stl');
    },
    downloadPLY: function() {
        if (ModelExporter.assertModelNonempty()) return;
        let ptr = Module.ccall('generatePLY', 'int', [], []);
        ModelExporter.downloadFile(ptr, ModelExporter.name+'.ply', 'model/stl');
    },
    downloadOBJ: function() {
        if (ModelExporter.assertModelNonempty()) return;
        let ptr = Module.ccall('generateOBJ', 'int', [], []);
        ModelExporter.downloadFile(ptr, ModelExporter.name+'.obj', 'model/obj');
    },
    downloadGLB: function() {
        if (ModelExporter.assertModelNonempty()) return;
        let ptr = Module.ccall('generateGLB', 'int', [], []);
        ModelExporter.downloadFile(ptr, ModelExporter.name+'.glb', 'model/gltf-binary');
    },
    init: function() {
        document.getElementById("export-stl").onclick = ModelExporter.downloadSTL;
        document.getElementById("export-ply").onclick = ModelExporter.downloadPLY;
        document.getElementById("export-obj").onclick = ModelExporter.downloadOBJ;
        document.getElementById("export-glb").onclick = ModelExporter.downloadGLB;
    }
};

function onError(message) {
    let errorMessage = document.getElementById("error-message");
    errorMessage.style.display = message ? "inline-block" : "none";
    errorMessage.innerHTML = message ? message : "";
}

function initDragDrop() {
    let dropArea = document.getElementById("control");
    let img = new Image();
    let onerror = function(message) {
        if (message)
            onError("error: " + message);
        else
            onError("error: unsupported image format");
    };
    onError();
    function imgOnload(img_src) {
        let canvas = document.createElement("canvas");
        if (img.width == 0 || img.height == 0)
            return onerror();
        var sc = Math.min(Math.sqrt(1048576/(img.width*img.height)), 1);
        canvas.width = img.width == 0 ? 1024 : Math.round(sc*img.width);
        canvas.height = img.height == 0 ? 1024 : Math.round(sc*img.height);
        let ctx = canvas.getContext("2d");
        ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
        let imageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
        let pixelData = imageData.data;
        let heapSpace = Module._malloc(pixelData.length);
        Module.HEAP8.set(pixelData, heapSpace);
        var name = img_src.split('/');
        name = name[name.length - 1].split('.');
        name = name.slice(0, name.length-1).join('.');
        ModelExporter.name = name;
        Module.ccall('updateImage',
            null, ['string', 'number', 'number', 'number'],
            [name, canvas.width, canvas.height, heapSpace]);
    }
    img.onload = function(e) { imgOnload(img.src); }
    img.src = "hermit_crab.svg";

    window.addEventListener("dragover", function(e) {
        e.preventDefault();
        dropArea.style.backgroundColor = "yellow";
    });
    window.addEventListener("dragleave", function(e) {
        e.preventDefault();
        dropArea.style.backgroundColor = null;
    });
    window.addEventListener("drop", function(e) {
        e.preventDefault();
        let errorMessage = document.getElementById("error-message");
        errorMessage.style.display = "none";
        dropArea.style.backgroundColor = null;
        if (e.dataTransfer.files.length > 0) {
            let file = e.dataTransfer.files[0];
            if (file && (
                file.type == "image/png" ||
                file.type == "image/jpeg" ||
                file.type == "image/gif" ||
                file.type == "image/webp" ||
                file.type == "image/svg+xml"
                )) {
                let reader = new FileReader();
                reader.onload = function(event) {
                    img = new Image();
                    img.onload = function(e) { imgOnload(file.name); };
                    img.onerror = onerror;
                    img.src = event.target.result;
                }
                reader.onerror = onerror;
                reader.readAsDataURL(file);
            }
            else onerror();
        }
        else if (e.dataTransfer.items.length > 0) {
            var found = false;
            for (var i = 0; i < e.dataTransfer.items.length; i++) {
                let item = e.dataTransfer.items[i]
                if (item.kind != "string" || item.type != "text/plain")
                    continue;
                item.getAsString((s) => {
                    let filename = s.split('/')[s.split('/').length-1];
                    img = new Image();
                    img.crossOrigin = 'anonymous';
                    img.onload = (e) => { imgOnload(filename); };
                    img.onerror = (e) => { onerror("failed to load image.<br/>" +
                        "If you are drag-dropping from browser, try switching browser, or save the image as a local file and drag drop from local file viewer."); };
                    img.src = s;
                });
                found = true;
            }
            if (!found) onerror();
        }
        else onerror();
    });
}

function calcScreenCenter() {
    let rect = document.getElementById("control").getBoundingClientRect();
    var w = window.innerWidth, h = window.innerHeight;
    var rl = rect.left, rb = h - rect.bottom;
    var cx = 0.5 * w, cy = 0.5 * h;
    if (rl > rb && rl > 0) cx = 0.5 * rl;
    else if (rb > 0) cy = 0.5 * rb;
    var com = { x: 2.0 * (cx / w - 0.5), y: 2.0 * (cy / h - 0.5) };
    com.x = Math.max(-0.6, Math.min(0.6, com.x));
    com.y = Math.max(-0.6, Math.min(0.6, com.y));
    return com;
}

function initInteraction() {

    // window resize
    let canvas = document.getElementById("emscripten-canvas");
    function onresize() {
        canvas.width = window.innerWidth;
        canvas.height = window.innerHeight;
        canvas.style.width = canvas.width + "px";
        canvas.style.height = canvas.height + "px";
        let com = calcScreenCenter();
        Module.ccall('resizeWindow',
            null,
            ['number', 'number', 'number', 'number'],
            [canvas.width, canvas.height,
                com.x, com.y]);
    };
    window.addEventListener("resize", onresize);
    onresize();

    // sliders - right click to reset
    let sliders = document.querySelectorAll('input[type=range]');
    for (var i = 0; i < sliders.length; i++) {
        sliders[i].addEventListener("contextmenu", (event) => {
            event.preventDefault();
            let slider = event.target;
            slider.value = Math.round((Number(slider.min)+Number(slider.max))/2);
            var event = new CustomEvent("input");
            slider.dispatchEvent(event);
        });
    }

}

function initConfig() {
    let selectChannel = document.getElementById("select-channel");
    let sliderThreshold = document.getElementById("slider-threshold");
    let checkboxReverse = document.getElementById("checkbox-reverse");
    let checkboxResampleBoundary = document.getElementById("checkbox-resample-boundary");
    let checkboxEdge = document.getElementById("checkbox-edge");
    let checkboxNormal = document.getElementById("checkbox-normal");
    let checkboxDoubleSided = document.getElementById("checkbox-double-sided");
    let checkboxTexture = document.getElementById("checkbox-texture");
    let sliderZPositive = document.getElementById("slider-z-pos");
    let sliderZNegative = document.getElementById("slider-z-neg");
    let checkboxClip = document.getElementById("checkbox-clip");
    let sliderClipZPositive = document.getElementById("slider-clip-z-pos");
    let sliderClipZNegative = document.getElementById("slider-clip-z-neg");
    function updateParameters() {
        onError();
        Module.ccall('setChannel', null, ['int'], [Number(selectChannel.value)]);
        Module.ccall('setThreshold', null, ['int'], [sliderThreshold.value]);
        Module.ccall('setAlphaReverse', null, ['int'], [checkboxReverse.checked]);
        Module.ccall('setResampleBoundary', null, ['int'], [checkboxResampleBoundary.checked]);
        Module.ccall('setMeshEdge', null, ['int'], [checkboxEdge.checked]);
        Module.ccall('setMeshNormal', null, ['int'], [checkboxNormal.checked]);
        Module.ccall('setMeshDoubleSided', null, ['int'], [checkboxDoubleSided.checked]);
        Module.ccall('setMeshTexture', null, ['int'], [checkboxTexture.checked]);
        Module.ccall('setZScale', null, ['number', 'number'], [
            Number(sliderZPositive.value) / (Number(sliderZPositive.max)+1),
            Number(sliderZNegative.value) / (Number(sliderZNegative.max)+1)
        ]);
        Module.ccall('setZClip', null, ['int', 'number', 'number'], [
            checkboxClip.checked,
            Number(sliderClipZPositive.value) / (Number(sliderClipZPositive.max)+1),
            Number(sliderClipZNegative.value) / (Number(sliderClipZNegative.max)+1)
        ]);
    }
    selectChannel.addEventListener("input", updateParameters);
    sliderThreshold.addEventListener("input", updateParameters);
    checkboxReverse.addEventListener("input", updateParameters);
    checkboxResampleBoundary.addEventListener("input", updateParameters);
    checkboxEdge.addEventListener("input", updateParameters);
    checkboxNormal.addEventListener("input", updateParameters);
    checkboxDoubleSided.addEventListener("input", updateParameters);
    checkboxTexture.addEventListener("input", updateParameters);
    sliderZPositive.addEventListener("input", updateParameters);
    sliderZNegative.addEventListener("input", updateParameters);
    checkboxClip.addEventListener("input", updateParameters);
    sliderClipZPositive.addEventListener("input", updateParameters);
    sliderClipZNegative.addEventListener("input", updateParameters);
    updateParameters();
    updateParameters();
}

function onReady() {
    initInteraction();
    initDragDrop();
    initConfig();
    ModelExporter.init();
}

