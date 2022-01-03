(function(){
    const SEC24H = 24*60*60.0;
    const DEG2RAD = Math.PI / 180.0;
    const SEC2RAD = 2*Math.PI / SEC24H;

    function clamp(x, lower, upper){
        return x < lower ? lower : x > upper ? upper : x;
    }

    //
    // Array Layout:
    // [0 4 8  12      [0
    //  1 5 9  13   *   1
    //  2 6 10 14       2
    //  3 7 11 15]      3]
    //
    const M4 = {
        identity: function(){
            return [1, 0, 0, 0,
                    0, 1, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1];
        },
        rotZ: function(rad){
            return [Math.cos(rad), Math.sin(rad), 0, 0,
                    -Math.sin(rad), Math.cos(rad), 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1];
        },
        rotY: function(rad){
            return [Math.cos(rad), 0, -Math.sin(rad), 0,
                    0, 1, 0, 0,
                    Math.sin(rad), 0, Math.cos(rad), 0,
                    0, 0, 0, 1];
        },
        rotX: function(rad){
            return [1, 0, 0, 0,
                    0, Math.cos(rad), Math.sin(rad), 0,
                    0, -Math.sin(rad),Math.cos(rad), 0,
                    0, 0, 0, 1];
        },
        perspective: function(fovYDeg, screenW, screenH, nearZ, farZ){
            const fovYRad = fovYDeg * DEG2RAD;
            const aspectRatio = screenH / screenW;

            const h = 1.0 / Math.tan(fovYRad / 2.0);
            const w = h * aspectRatio;
            const zNearFar = nearZ - farZ;
            return [
                w, 0, 0, 0,
                0, h, 0, 0,
                0, 0, (farZ+nearZ)/zNearFar, -1,
                0, 0, 2*nearZ*farZ/zNearFar, 0];
        },
        translate: function(dx, dy, dz){
            return [
                1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                dx, dy, dz, 1];
        },
        scale: function(sx, sy, sz){
            if(sy === undefined){sy = sx;}
            if(sz === undefined){sz = sx;}
            return [
                sx, 0, 0, 0,
                0, sy, 0, 0,
                0, 0, sz, 0,
                0, 0, 0, 1];
        },
        mulM4: function(l, r){
            return [
                l[0]*r[0]+l[4]*r[1]+l[8]*r[2]+l[12]*r[3],
                l[1]*r[0]+l[5]*r[1]+l[9]*r[2]+l[13]*r[3],
                l[2]*r[0]+l[6]*r[1]+l[10]*r[2]+l[14]*r[3],
                l[3]*r[0]+l[7]*r[1]+l[11]*r[2]+l[15]*r[3],

                l[0]*r[4]+l[4]*r[5]+l[8]*r[6]+l[12]*r[7],
                l[1]*r[4]+l[5]*r[5]+l[9]*r[6]+l[13]*r[7],
                l[2]*r[4]+l[6]*r[5]+l[10]*r[6]+l[14]*r[7],
                l[3]*r[4]+l[7]*r[5]+l[11]*r[6]+l[15]*r[7],

                l[0]*r[8]+l[4]*r[9]+l[8]*r[10]+l[12]*r[11],
                l[1]*r[8]+l[5]*r[9]+l[9]*r[10]+l[13]*r[11],
                l[2]*r[8]+l[6]*r[9]+l[10]*r[10]+l[14]*r[11],
                l[3]*r[8]+l[7]*r[9]+l[11]*r[10]+l[15]*r[11],

                l[0]*r[12]+l[4]*r[13]+l[8]*r[14]+l[12]*r[15],
                l[1]*r[12]+l[5]*r[13]+l[9]*r[14]+l[13]*r[15],
                l[2]*r[12]+l[6]*r[13]+l[10]*r[14]+l[14]*r[15],
                l[3]*r[12]+l[7]*r[13]+l[11]*r[14]+l[15]*r[15]];
        },
        mulV4: function(l, r){
            return [
                l[0] * r[0] + l[4] * r[1] + l[8] * r[2] + l[12] * r[3],
                l[1] * r[0] + l[5] * r[1] + l[9] * r[2] + l[13] * r[3],
                l[2] * r[0] + l[6] * r[1] + l[10] * r[2] + l[14] * r[3],
                l[3] * r[0] + l[7] * r[1] + l[11] * r[2] + l[15] * r[3]];
        }
    };

    function dirRADec(ra, dec){
        // x+:6h, y+: North, z+:0h(Vernal Equinox)
        return [
            Math.cos(dec) * Math.sin(ra),
            Math.sin(dec),
            Math.cos(dec) * Math.cos(ra),
            1];
    }
    function dirAzEl(az, el){
        // x+:East, y+:Zenith, z+:South
        return [
            Math.cos(el)*Math.sin(az),
            Math.sin(el),
            Math.cos(el)*Math.cos(az),
            1];
    }

    function convertRADecToAzEl(ra, dec, matEquToHor)
    {
        const dirEquatorial = dirRADec(ra, dec);
        const dirHorizontal = M4.mulV4(matEquToHor, dirEquatorial);
        const x = dirHorizontal[0];//for East
        const y = dirHorizontal[1];//for Zenith
        const z = dirHorizontal[2];//for South
        const az = Math.atan2(-x, z);
        const el = Math.asin(y);

        return {
            az: az * 180 / Math.PI,
            el: el * 180 / Math.PI
        };
    }

    //
    // Time Conversion
    //
    function convertUnixSecondsToGreenwichMeanSiderealTime(unixSeconds){
        ///@todo validation
        ///@todo add leap seconds
        const jd00 = unixSeconds / SEC24H + (2440587.5 - 2451545.0); //2440587.5=Unix Epoch(in JD), 2451545.0=J2000.0(in JD)
        const t = jd00 / 36525.0; //36525.0=Days per Julian century
        const f = SEC24H * (jd00 % 1.0);
        const A = 24110.54841  -  SEC24H / 2.0;
        const B = 8640184.812866;
        const C = 0.093104;
        const D =  -6.2e-6;
        const gmst = ((A + (B + (C + D * t) * t) * t) + f) * SEC2RAD; //[rad]
        const gmstNormalized = gmst % (2*Math.PI);
        return gmstNormalized < 0 ? (2*Math.PI) + gmstNormalized : gmstNormalized;
    }

    //
    // StarChart
    //

    function StarChart(options){
        // astronomical latitude and longitude
        const lat = typeof(options.lat)=="number" ? options.lat : 35.681236; //Tokyo Station(geodetic)
        const lng = typeof(options.lng)=="number" ? options.lng : 139.767125; //Tokyo Station(geodetic)

        var matEquToHor; //x+:6h, y+:North, z+:0h to x+:East, y+:Zenith, z+:South
        function updateEquToHorMatrix()
        {
            const st = convertUnixSecondsToGreenwichMeanSiderealTime(Math.floor(Date.now() / 1000)) + lng*DEG2RAD;
            matEquToHor = M4.mulM4(M4.rotX((lat-90) * DEG2RAD), M4.rotY(-st));
        }
        updateEquToHorMatrix();

        //
        // View
        //

        let screenWidth;
        let screenHeight;
        let fovY = typeof(options.fov)=="number" ? options.fov : 50; //[deg]
        const useStereographicProjection = true;
        var matProj;
        function updateProjectionMatrix()
        {
            matProj = M4.perspective(fovY, screenWidth, screenHeight, 0.0125, 3.0);
            if(useStereographicProjection){
                matProj = M4.mulM4(matProj, M4.translate(0, 0, -1)); //Stereographic Projection
            }
            updateViewMatrix();
        }
        function setFOV(fov)
        {
            fov = clamp(fov, 30, 120);
            if(fov != fovY){
                fovY = fov;
                updateProjectionMatrix();
                drawFrame();
            }
        }
        function setFOVDelta(delta)
        {
            setFOV(fovY + delta);
        }

        if(typeof(options.ra)=="number" && typeof(options.dec)=="number"){
            const horDir = convertRADecToAzEl(options.ra*DEG2RAD, options.dec*DEG2RAD, matEquToHor);
            options.az = horDir.az;
            options.el = horDir.el;
        }
        let viewAz = typeof(options.az)=="number" ? options.az : 0;
        let viewEl = typeof(options.el)=="number" ? options.el : 0;
        let matView;
        function updateViewMatrix(){
            matView = M4.mulM4(M4.rotX(-viewEl * DEG2RAD), M4.rotY((180+viewAz) * DEG2RAD));
        }
        function setViewDir(az, el){
            viewAz = az;
            viewEl = clamp(el, -120, 120); //more than 90 degrees to see behind
            updateViewMatrix();
            drawFrame();
        }

        let cv;
        let renderer;
        function setupCanvas()
        {
            cv = document.createElement("canvas");
            cv.style = "position: absolute; left:0; top:0;";
            document.body.appendChild(cv);

            window.addEventListener("resize", onResizeWindow);
            function onResizeWindow(){
                cv.width = screenWidth = window.innerWidth;
                cv.height = screenHeight = window.innerHeight;
                updateProjectionMatrix();
                if(renderer){
                    renderer.onResize();
                    drawFrame();
                }
            }
            onResizeWindow();

            renderer = new WebGLRenderer(cv, onRendererIsReady);
            function onRendererIsReady(){
                drawFrame();
            }
        }

        function drawFrame()
        {
            if(renderer && renderer.isReady()){
                renderer.drawFrame(matEquToHor, matView, matProj);
            }
        }


        // setup
        updateEquToHorMatrix();
        setupCanvas();

        // public
        this.getCV = function(){return cv;};
        this.getFOV = function(){return fovY;};
        this.setFOV = setFOV;
        this.zoom = setFOVDelta;
        this.getViewAz = function(){return viewAz;};
        this.getViewEl = function(){return viewEl;};
        this.setViewDir = setViewDir;
    }//StarChart

    //
    // Renderer
    //

    const COLOR_GRID = [0, 0.5, 1, 0.25];
    const COLOR_HORIZON = [0.5, 1, 0, 0.4];
    const COLOR_NORTH = [1, 0, 0, 0.4];
    const COLOR_SOUTH = [1, 1, 1, 0.25];

    function WebGLRenderer(cv, callbackReady)
    {
        const gl = cv.getContext("webgl");

        this.drawFrame = drawFrame;
        function drawFrame(matEquToHor, matView, matProj)
        {
            clear();
            drawCelestialSphere(matEquToHor, matView, matProj);
            drawGrid(matEquToHor, matView, matProj);
        }

        function clear()
        {
            gl.clearColor(0, 0, 0, 1);
            gl.clear(gl.COLOR_BUFFER_BIT);
        }


        function ArcRenderer(){
            const step = 5;
            const numVertices = 360 / step + 1;
            const bufferLayout = {
                numVertices: numVertices,
                attrStride: 2 * Float32Array.BYTES_PER_ELEMENT,
                attr: {
                    position: {
                        offset: 0,
                        size: 2,
                    },
                }
            };
            const bufferData = new Float32Array(numVertices * 2);
            for(let i = 0; i < numVertices; ++i){
                const rad = i*step*DEG2RAD;
                bufferData[i*2+0] = Math.cos(rad);
                bufferData[i*2+1] = Math.sin(rad);
            }
            const buffer = gl.createBuffer();
            gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
            gl.bufferData(gl.ARRAY_BUFFER, bufferData, gl.STATIC_DRAW);

            var shader = createShaderProgram(
                "attribute vec2 a_position; uniform mat4 u_matrix;"+
                    "void main(void){gl_Position = u_matrix * vec4(a_position, 0.0, 1.0);}",
                "precision mediump float; uniform vec4 u_color; void main(void){gl_FragColor = u_color;}");
            const shaderLocations = {
                attr: {
                    position: gl.getAttribLocation(shader, "a_position")
                },
                uniform: {
                    matrix: gl.getUniformLocation(shader, "u_matrix"),
                    color: gl.getUniformLocation(shader, "u_color")
                },
            };

            function draw(mat, color, angleBegin, angleEnd){
                if(angleBegin === undefined){ angleBegin = 0;}
                if(angleEnd === undefined){ angleEnd = 360;}
                angleBegin = clamp(Math.floor(angleBegin / step), 0, numVertices);
                angleEnd = clamp(Math.floor(angleEnd / step) + 1, 0, numVertices);

                gl.useProgram(shader);
                gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
                gl.vertexAttribPointer(
                    shaderLocations.attr.position,
                    bufferLayout.attr.position.size,
                    gl.FLOAT, false, bufferLayout.attrStride,
                    bufferLayout.attr.position.offset);
                gl.enableVertexAttribArray(shaderLocations.attr.position);
                gl.uniformMatrix4fv(shaderLocations.uniform.matrix, false, new Float32Array(mat));
                gl.uniform4fv(shaderLocations.uniform.color, new Float32Array(color));
                gl.drawArrays(gl.LINE_STRIP, angleBegin, angleEnd - angleBegin);
                gl.disableVertexAttribArray(shaderLocations.attr.position);
            }
            this.draw = draw;
        }
        const arcRenderer = new ArcRenderer();
        function drawGrid(matEquToHor, matView, matProj){
            gl.enable(gl.BLEND);
            gl.blendFuncSeparate(gl.SRC_ALPHA, gl.ONE_MINUS_SOURCE_ALPHA, gl.ONE, gl.ONE);

            const matHorToProj = M4.mulM4(matProj, matView);
            for(let el = -90+15; el < 90; el += 15){
                arcRenderer.draw(
                    M4.mulM4(matHorToProj, M4.mulM4(M4.mulM4(M4.translate(0, Math.sin(el*DEG2RAD), 0), M4.scale(Math.cos(el*DEG2RAD))), M4.rotX(Math.PI/2))),
                    el==0 ? COLOR_HORIZON : COLOR_GRID);
            }
            for(let az = 0; az < 180; az += 15){
                const mat = M4.mulM4(matHorToProj, M4.mulM4(M4.rotY(az*DEG2RAD), M4.rotZ(90*DEG2RAD)));
                if(az == 90){
                    arcRenderer.draw(mat, COLOR_SOUTH, 0, 180);
                    arcRenderer.draw(mat, COLOR_NORTH, 180, 360);
                }
                else{
                    arcRenderer.draw(mat, COLOR_GRID);
                }
            }
            gl.disable(gl.BLEND);
        }

        function CelestialSphereRenderer(filename, callbackOnLoad)
        {
            const DIV_SPHERE = 32;
            const vertices = new Float32Array((DIV_SPHERE+1) * (DIV_SPHERE * 2 + 1) * 5);
            for(let i = 0, vi = 0; i <= DIV_SPHERE; ++i){
                const ai = i * (Math.PI / DIV_SPHERE);
                // x+:6h, y+: North, z+:0h(Vernal Equinox)
                const y = Math.cos(ai);
                const r = Math.sin(ai);
                for(let j = 0; j <= 2*DIV_SPHERE; ++j){
                    const aj = j * (Math.PI / DIV_SPHERE);
                    const x = r * Math.sin(aj);
                    const z = r * Math.cos(aj);
                    // Left:24h - Right:0h
                    // Top:North - Bottom:South
                    const s = 1.0 - j / (2*DIV_SPHERE);
                    const t = i / DIV_SPHERE;
                    vertices[vi++] = x;
                    vertices[vi++] = y;
                    vertices[vi++] = z;
                    vertices[vi++] = s;
                    vertices[vi++] = t;
                }
            }
            const numIndices = 1 + DIV_SPHERE * (1 + 2*DIV_SPHERE * 2);
            const indices = new Uint16Array(numIndices);
            let ii = 0;
            let dir = 1;
            indices[ii++] = 0;
            for(let i = 0; i < DIV_SPHERE; ++i){
                const indexTop = i * (2*DIV_SPHERE + 1);
                const indexBottom = indexTop + (2*DIV_SPHERE + 1);
                const indexBottomRight = indexBottom + (dir >= 0 ? 0 : 2*DIV_SPHERE);
                indices[ii++] = indexBottomRight;
                for(let j = 1; j <= 2*DIV_SPHERE; ++j){
                    const indexTopLeft = indexTop + (dir >= 0 ? j : (2*DIV_SPHERE - j));
                    const indexBottomLeft = indexBottom + (dir >= 0 ? j : (2*DIV_SPHERE - j));
                    indices[ii++] = indexTopLeft;
                    indices[ii++] = indexBottomLeft;
                }
                dir = -dir;
            }
            if(ii != numIndices){
                throw Error("Assertion failed (ii != numIndices)");
            }
            const vertexBuffer = gl.createBuffer();
            gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
            gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);
            const indexBuffer = gl.createBuffer();
            gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);
            gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, indices, gl.STATIC_DRAW);
            const bufferLayout = {
                numVertices: (DIV_SPHERE+1) * (DIV_SPHERE*2+1),
                attrStride: 5 * Float32Array.BYTES_PER_ELEMENT,
                attr: {
                    position: {
                        offset: 0 * Float32Array.BYTES_PER_ELEMENT,
                        size: 3,
                    },
                    texcoord: {
                        offset: 3 * Float32Array.BYTES_PER_ELEMENT,
                        size: 2,
                    },
                }
            };

            const shader = createShaderProgram(
                "attribute vec3 a_position; attribute vec2 a_texcoord; uniform mat4 u_matrix; varying highp vec2 v_texcoord;"+
                    "void main(void){gl_Position = u_matrix * vec4(a_position, 1.0); v_texcoord = a_texcoord;}",
                "varying highp vec2 v_texcoord; uniform sampler2D u_sampler;" +
                    "void main(void){gl_FragColor = texture2D(u_sampler, vec2(v_texcoord.s, v_texcoord.t));}");
            const shaderLocations = {
                attr: {
                    position: gl.getAttribLocation(shader, "a_position"),
                    texcoord: gl.getAttribLocation(shader, "a_texcoord"),
                },
                uniform: {
                    matrix: gl.getUniformLocation(shader, "u_matrix"),
                    sampler: gl.getUniformLocation(shader, "u_sampler")
                },
            };

            let texture;
            const texImage = new Image();
            texImage.onload = onLoadTexture;
            texImage.src = filename;
            function onLoadTexture(){
                texture = gl.createTexture();
                gl.bindTexture(gl.TEXTURE_2D, texture);
                gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texImage);
                gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
                const useMipMap = false;
                if(useMipMap){
                    gl.generateMipmap(gl.TEXTURE_2D);
                    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MAPMAP_NEAREST);
                    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT);
                    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
                }
                else{
                    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
                    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
                    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
                }
                gl.bindTexture(gl.TEXTURE_2D, null);

                if(callbackOnLoad){
                    callbackOnLoad();
                }
            }

            function draw(mat){
                if(!shader || !texture){
                    return;
                }
                gl.useProgram(shader);
                gl.bindBuffer(gl.ARRAY_BUFFER, vertexBuffer);
                gl.vertexAttribPointer(
                    shaderLocations.attr.position,
                    bufferLayout.attr.position.size,
                    gl.FLOAT, false, bufferLayout.attrStride,
                    bufferLayout.attr.position.offset);
                gl.vertexAttribPointer(
                    shaderLocations.attr.texcoord,
                    bufferLayout.attr.texcoord.size,
                    gl.FLOAT, false, bufferLayout.attrStride,
                    bufferLayout.attr.texcoord.offset);
                gl.enableVertexAttribArray(shaderLocations.attr.position);
                gl.enableVertexAttribArray(shaderLocations.attr.texcoord);
                gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, indexBuffer);
                gl.uniformMatrix4fv(shaderLocations.uniform.matrix, false, new Float32Array(mat));
                gl.activeTexture(gl.TEXTURE0);
                gl.bindTexture(gl.TEXTURE_2D, texture);
                gl.uniform1i(shaderLocations.uniform.sampler, 0);
                gl.drawElements(gl.TRIANGLE_STRIP, numIndices, gl.UNSIGNED_SHORT, 0);
                gl.disableVertexAttribArray(shaderLocations.attr.position);
                gl.disableVertexAttribArray(shaderLocations.attr.texcoord);
                gl.bindTexture(gl.TEXTURE_2D, null);
                gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
                gl.bindBuffer(gl.ARRAY_BUFFER, null);
            }
            this.draw = draw;
        }
        const celestialSphereRenderer = new CelestialSphereRenderer("equirectangular.jpg", onCelestialSphereLoaded);
        function drawCelestialSphere(matEquToHor, matView, matProj)
        {
            const mat = M4.mulM4(M4.mulM4(matProj, matView), matEquToHor);
            celestialSphereRenderer.draw(mat);
        }
        let celestialSphereLoaded = false;
        function onCelestialSphereLoaded()
        {
            celestialSphereLoaded = true;
            if(callbackReady){
                callbackReady();
            }
        }
        function isReady()
        {
            return celestialSphereLoaded;
        }
        this.isReady = isReady;

        //
        // WebGL Util
        //
        function compileShader(source, type)
        {
            var shader = gl.createShader(type); //gl.FRAGMENT_SHADER, VERTEX_SHADER
            gl.shaderSource(shader, source);
            gl.compileShader(shader);
            if(!gl.getShaderParameter(shader, gl.COMPILE_STATUS)){
                throw new Error("Failed to compile shader." + gl.getShaderInfoLog(shader));
            }
            return shader;
        }
        function createShaderProgram(vertexShaderSource, fragmentShaderSource)
        {
            var fragmentShader = compileShader(fragmentShaderSource, gl.FRAGMENT_SHADER);
            var vertexShader = compileShader(vertexShaderSource, gl.VERTEX_SHADER);
            var program = gl.createProgram();
            gl.attachShader(program, vertexShader);
            gl.attachShader(program, fragmentShader);
            gl.linkProgram(program);
            if(!gl.getProgramParameter(program, gl.LINK_STATUS)){
                throw new Error("Failed to link shader program. ");
            }
            return program;
        }

        function onResize(){
            gl.viewport(0, 0, cv.width, cv.height);
        }
        this.onResize = onResize;
    }//WebGLRenderer




    //
    // Control
    //

    function StarChartController(starChart){
        const cv = starChart.getCV();

        //
        // Touch
        //
        const currentTouches = [];
        function findCurrentTouchById(id){
            return currentTouches.findIndex((t)=>t.id==id);
        }
        function averageXY(touches){
            let x = 0, y = 0;
            for(let i = 0; i < touches.length; ++i){
                x += touches[i].x;
                y += touches[i].y;
            }
            return {x: x / touches.length,
                    y: y / touches.length};
        }

        let touchStartState;
        function refreshTouchStartState(){
            touchStartState = {
                fov: starChart.getFOV(),
                viewAz: starChart.getViewAz(),
                viewEl: starChart.getViewEl(),
                touches: currentTouches.slice(),
                center: averageXY(currentTouches)
            };
        }

        cv.addEventListener("touchstart", onTouchStart, false);
        cv.addEventListener("touchend", onTouchEnd, false);
        cv.addEventListener("touchcancel", onTouchCancel, false);
        cv.addEventListener("touchmove", onTouchMove, false);
        function onTouchStart(ev){
            ev.preventDefault();
            const touches = ev.changedTouches;
            for(let ti = 0; ti < touches.length; ++ti){
                currentTouches.push({
                    id: touches[ti].identifier,
                    x: touches[ti].clientX,
                    y: touches[ti].clientY,
                });
            }
            refreshTouchStartState();
        }
        function onTouchEnd(ev){
            ev.preventDefault();
            const touches = ev.changedTouches;
            for(let ti = 0; ti < touches.length; ++ti){
                const index = findCurrentTouchById(touches[ti].identifier);
                if(index >= 0){
                    currentTouches.splice(index);
                }
            }
            refreshTouchStartState();
        }
        function onTouchCancel(ev){
            onTouchEnd(ev);
        }
        function onTouchMove(ev){
            ev.preventDefault();
            const touches = ev.changedTouches;
            for(let ti = 0; ti < touches.length; ++ti){
                const index = findCurrentTouchById(touches[ti].identifier);
                if(index >= 0){
                    currentTouches[index] = {
                        id: touches[ti].identifier,
                        x: touches[ti].clientX,
                        y: touches[ti].clientY,
                    };
                }
            }

            const anglePerPixel = 120 / cv.height;
            // Move
            if(currentTouches.length >= 1 && touchStartState.touches.length >= 1){
                const currCenter = averageXY(currentTouches);
                const dx = (currCenter.x - touchStartState.center.x) * anglePerPixel;
                const dy = (currCenter.y - touchStartState.center.y) * anglePerPixel;
                starChart.setViewDir(
                    touchStartState.viewAz - dx,
                    touchStartState.viewEl + dy);
            }

            // Zoom
            if(currentTouches.length >= 2 && touchStartState.touches.length >= 2){
                const startDistanceX = touchStartState.touches[0].x - touchStartState.touches[1].x;
                const startDistanceY = touchStartState.touches[0].y - touchStartState.touches[1].y;
                const startDistance = Math.sqrt(startDistanceX*startDistanceX + startDistanceY*startDistanceY);
                const currDistanceX = currentTouches[0].x - currentTouches[1].x;
                const currDistanceY = currentTouches[0].y - currentTouches[1].y;
                const currDistance = Math.sqrt(currDistanceX*currDistanceX + currDistanceY*currDistanceY);
                const delta = (currDistance - startDistance) * anglePerPixel;
                starChart.setFOV(touchStartState.fov - delta);
            }
        }

        // Mouse
        cv.addEventListener("mousedown", onMouseDown, false);
        cv.addEventListener("mouseup", onMouseUp, false);
        cv.addEventListener("mousemove", onMouseMove, false);
        cv.addEventListener("wheel", onMouseWheel, false);
        var mouseDownState = null;
        function onMouseDown(ev){
            ev.preventDefault();
            mouseDownState = {x:ev.clientX, y:ev.clientY, viewAz:starChart.getViewAz(), viewEl:starChart.getViewEl()};
        }
        function onMouseUp(ev){
            ev.preventDefault();
            mouseDownState = null;
        }
        function onMouseMove(ev){
            const anglePerPixel = 2 * starChart.getFOV() / cv.height;
            if(mouseDownState){
                ev.preventDefault();
                starChart.setViewDir(
                    mouseDownState.viewAz - (ev.clientX - mouseDownState.x) * anglePerPixel,
                    mouseDownState.viewEl + (ev.clientY - mouseDownState.y) * anglePerPixel);
            }
        }
        function onMouseWheel(ev){
            ev.preventDefault();
            starChart.zoom(ev.deltaY * 0.01);
        }
    }


    function setup(options){
        const starChart = new StarChart(options);
        const controller = new StarChartController(starChart);
    }

    function parseQueryString()
    {
        const options = {};
        const q = document.location.search.substr(1);
        if(q.length > 0){
            const ps = q.split("&");
            for(let pi = 0; pi < ps.length; ++pi){
                const kv = ps[pi].split("=");
                const key = kv[0];
                const value = decodeURI(kv[1].replace(/\+/g, " "));
                switch(key){
                case "fov":
                case "ra":
                case "dec":
                case "az":
                case "el":
                case "lng":
                case "lat":
                    options[key] = parseFloat(value); break;
                case "renderer":
                    options[key] = value; break;
                }
            }
        }
        return options;
    }

    function main(){
        setup(parseQueryString());
    }

    window.SkyView = {
        main: main
    };
})();
