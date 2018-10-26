
(function(d, f) {
    var a = d.jQuery || d.Cowboy || (d.Cowboy = {}), e;
    a.throttle = e = function(i, c, b, g) {
        function d() {
            function a() {
                e = +new Date;
                b.apply(k, l)
            }
            function m() {
                h = f
            }
            var k = this, j = +new Date - e, l = arguments;
            g && !h && a();
            h && clearTimeout(h);
            g === f && j > i ? a() : c !== true && (h = setTimeout(g ? m : a, g === f ? i - j : i))
        }
        var h, e = 0;
        typeof c !== "boolean" && (g = b, b = c, c = f);
        if (a.guid)
            d.guid = b.guid = b.guid || a.guid++;
        return d
    };
    a.debounce = function(a, c, b) {
        return b === f ? e(a, c, false) : e(a, b, c !== false)
    }
})(this);
(function($) {
    $.idleTimer = function(newTimeout, elem) {
        var idle = false, enabled = true, timeout = 30000, events = 'mousemove keydown DOMMouseScroll mousewheel mousedown touchstart touchmove';
        elem = elem || document;
        var toggleIdleState = function(myelem) {
            if (typeof myelem === 'number') {
                myelem = undefined;
            }
            var obj = $.data(myelem || elem, 'idleTimerObj');
            obj.idle = !obj.idle;
            var elapsed = (+new Date()) - obj.olddate;
            obj.olddate = +new Date();
            if (obj.idle && (elapsed < timeout)) {
                obj.idle = false;
                clearTimeout($.idleTimer.tId);
                if (enabled)
                    $.idleTimer.tId = setTimeout(toggleIdleState, timeout);
                return;
            }
            var event = jQuery.Event($.data(elem, 'idleTimer', obj.idle ? "idle" : "active") + '.idleTimer');
            $(elem).trigger(event);
        }, stop = function(elem) {
            var obj = $.data(elem, 'idleTimerObj');
            obj.enabled = false;
            clearTimeout(obj.tId);
            $(elem).off('.idleTimer');
        }, handleUserEvent = function() {
            var obj = $.data(this, 'idleTimerObj');
            clearTimeout(obj.tId);
            if (obj.enabled) {
                if (obj.idle) {
                    toggleIdleState(this);
                }
                obj.tId = setTimeout(toggleIdleState, obj.timeout);
            }
        };
        var obj = $.data(elem, 'idleTimerObj') || {};
        obj.olddate = obj.olddate || +new Date();
        if (typeof newTimeout === "number") {
            timeout = newTimeout;
        } else if (newTimeout === 'destroy') {
            stop(elem);
            return this;
        } else if (newTimeout === 'getElapsedTime') {
            return(+new Date()) - obj.olddate;
        }
        $(elem).on($.trim((events + ' ').split(' ').join('.idleTimer ')), handleUserEvent);
        obj.idle = idle;
        obj.enabled = enabled;
        obj.timeout = timeout;
        obj.tId = setTimeout(toggleIdleState, obj.timeout);
        $.data(elem, 'idleTimer', "active");
        $.data(elem, 'idleTimerObj', obj);
    };
    $.fn.idleTimer = function(newTimeout) {
        if (this[0]) {
            $.idleTimer(newTimeout, this[0]);
        }
        return this;
    };
})(jQuery);
(function(a, b, c) {
    function z(a, b, c) {
        var e = jQuery.inArray(a, d) > -1;
        (a == "width" || a == "height") && b === parseFloat(c.css(a)) && (e = !1);
        return e
    }
    function y(a) {
        return parseFloat(a.replace(/px/i, ""))
    }
    function x(a) {
        for (var b in a)
            return!1;
        return!0
    }
    function w(a) {
        for (var b in a)
            if ((b == "width" || b == "height") && (a[b] == "show" || a[b] == "hide" || a[b] == "toggle"))
                return!0;
        return!1
    }
    function v(a, b, c, d, e, g, h, i) {
        a = typeof a == "undefined" ? {} : a, a.secondary = typeof a.secondary == "undefined" ? {} : a.secondary;
        for (var j = f.length - 1; j >= 0; j--)
            typeof a[f[j] + "transition-property"] == "undefined" && (a[f[j] + "transition-property"] = ""), a[f[j] + "transition-property"] += ", " + (g === !0 && h === !0 ? f[j] + "transform" : b), a[f[j] + "transition-duration"] = c + "ms", a[f[j] + "transition-timing-function"] = d, a.secondary[g === !0 && h === !0 ? f[j] + "transform" : b] = g === !0 && h === !0 ? t(a.meta.left, a.meta.top, i) : e;
        return a
    }
    function u(a, b, c, d, f, g, h, i) {
        var l = a.data(k) ? x(a.data(k)) ? jQuery.extend(!0, {}, j) : a.data(k) : jQuery.extend(!0, {}, j), m = f, n = jQuery.inArray(b, e) > -1;
        if (n) {
            var o = l.meta, p = y(a.css(b)) || 0, q = b + "_o";
            m = n ? f - p : f, o[b] = m, o[q] = a.css(b) == "auto" ? 0 + m : p + m || 0, l.meta = o, h && m === 0 && (m = 0 - o[q], o[b] = m, o[q] = 0)
        }
        return a.data(k, v(l, b, c, d, m, g, h, i))
    }
    function t(a, b, c) {
        return c === !0 && r ? "translate3d(" + a + "px," + b + "px,0)" : "translate(" + a + "px," + b + "px)"
    }
    function s(a, b, c, d) {
        var e = h.exec(b), f = a.css(c) === "auto" ? 0 : a.css(c), g = typeof f == "string" ? y(f) : f, i = typeof b == "string" ? y(b) : b, j = d === !0 ? 0 : g, k = a.is(":hidden"), l = a.translation();
        c == "left" && (j = parseInt(g, 10) + l.x), c == "top" && (j = parseInt(g, 10) + l.y), !e && b == "show" ? (j = 1, k && a.css({display: "block", opacity: 0})) : !e && b == "hide" && (j = 0);
        if (e) {
            var m = parseFloat(e[2]);
            e[1] && (m = (e[1] === "-=" ? -1 : 1) * m + parseInt(j, 10));
            return m
        }
        return j
    }
    var d = ["top", "right", "bottom", "left", "opacity", "height", "width"], e = ["top", "right", "bottom", "left"], f = ["", "-webkit-", "-moz-", "-o-"], g = ["avoidTransforms", "useTranslate3d", "leaveTransforms"], h = /^([+-]=)?([\d+-.]+)(.*)$/, i = /([A-Z])/g, j = {secondary: {}, meta: {top: 0, right: 0, bottom: 0, left: 0}}, k = "jQe", l = "cubic-bezier(", m = ")", n = document.body || document.documentElement, o = n.style, p = o.WebkitTransition !== undefined ? "webkitTransitionEnd" : o.OTransition !== undefined ? "oTransitionEnd" : "transitionend", q = o.WebkitTransition !== undefined || o.MozTransition !== undefined || o.OTransition !== undefined || o.transition !== undefined, r = "WebKitCSSMatrix"in window && "m11"in new WebKitCSSMatrix;
    jQuery.fn.translation = function() {
        if (!this[0])
            return null;
        var a = this[0], b = window.getComputedStyle(a, null), c = {x: 0, y: 0};
        for (var d = f.length - 1; d >= 0; d--) {
            var e = b.getPropertyValue(f[d] + "transform");
            if (e && /matrix/i.test(e)) {
                var g = e.replace(/^matrix\(/i, "").split(/, |\)$/g);
                c = {x: parseInt(g[4], 10), y: parseInt(g[5], 10)};
                break
            }
        }
        return c
    }, jQuery.fn.animate = function(a, c, d, h) {
        a = a || {};
        var i = typeof a.bottom == "undefined" && typeof a.right == "undefined", j = jQuery.speed(c, d, h), n = this, o = 0, r = function() {
            o--, o === 0 && typeof j.complete == "function" && j.complete.apply(n[0], arguments)
        };
        if (!q || x(a) || w(a) || j.duration <= 0 || jQuery.fn.animate.defaults.avoidTransforms === !0 && a.avoidTransforms !== !1)
            return b.apply(this, arguments);
        return this[j.queue === !0 ? "queue" : "each"](function() {
            var c = jQuery(this), d = jQuery.extend({}, j), h = function() {
                var b = {};
                for (var d = f.length - 1; d >= 0; d--)
                    b[f[d] + "transition-property"] = "none", b[f[d] + "transition-duration"] = "", b[f[d] + "transition-timing-function"] = "";
                c.unbind(p);
                if (!a.leaveTransforms == !0) {
                    var g = c.data(k) || {}, h = {};
                    for (d = f.length - 1; d >= 0; d--)
                        h[f[d] + "transform"] = "";
                    if (i && typeof g.meta != "undefined")
                        for (var j = 0, l; l = e[j]; ++j)
                            h[l] = g.meta[l + "_o"] + "px";
                    c.css(b).css(h)
                }
                a.opacity === "hide" && c.css("display", "none"), c.data(k, null), r.call(c)
            }, n = {bounce: l + "0.0, 0.35, .5, 1.3" + m, linear: "linear", swing: "ease-in-out", easeInQuad: l + "0.550, 0.085, 0.680, 0.530" + m, easeInCubic: l + "0.550, 0.055, 0.675, 0.190" + m, easeInQuart: l + "0.895, 0.030, 0.685, 0.220" + m, easeInQuint: l + "0.755, 0.050, 0.855, 0.060" + m, easeInSine: l + "0.470, 0.000, 0.745, 0.715" + m, easeInExpo: l + "0.950, 0.050, 0.795, 0.035" + m, easeInCirc: l + "0.600, 0.040, 0.980, 0.335" + m, easeOutQuad: l + "0.250, 0.460, 0.450, 0.940" + m, easeOutCubic: l + "0.215, 0.610, 0.355, 1.000" + m, easeOutQuart: l + "0.165, 0.840, 0.440, 1.000" + m, easeOutQuint: l + "0.230, 1.000, 0.320, 1.000" + m, easeOutSine: l + "0.390, 0.575, 0.565, 1.000" + m, easeOutExpo: l + "0.190, 1.000, 0.220, 1.000" + m, easeOutCirc: l + "0.075, 0.820, 0.165, 1.000" + m, easeInOutQuad: l + "0.455, 0.030, 0.515, 0.955" + m, easeInOutCubic: l + "0.645, 0.045, 0.355, 1.000" + m, easeInOutQuart: l + "0.770, 0.000, 0.175, 1.000" + m, easeInOutQuint: l + "0.860, 0.000, 0.070, 1.000" + m, easeInOutSine: l + "0.445, 0.050, 0.550, 0.950" + m, easeInOutExpo: l + "1.000, 0.000, 0.000, 1.000" + m, easeInOutCirc: l + "0.785, 0.135, 0.150, 0.860" + m}, q = {}, t = n[d.easing || "swing"] ? n[d.easing || "swing"] : d.easing || "swing";
            for (var v in a)
                if (jQuery.inArray(v, g) === -1) {
                    var w = jQuery.inArray(v, e) > -1, y = s(c, a[v], v, w && a.avoidTransforms !== !0);
                    a.avoidTransforms !== !0 && z(v, y, c) ? u(c, v, d.duration, t, w && a.avoidTransforms === !0 ? y + "px" : y, w && a.avoidTransforms !== !0, i, a.useTranslate3d === !0) : q[v] = a[v]
                }
            var A = c.data(k) || {};
            for (var B = f.length - 1; B >= 0; B--)
                typeof A[f[B] + "transition-property"] != "undefined" && (A[f[B] + "transition-property"] = A[f[B] + "transition-property"].substr(2));
            c.data(k, A).unbind(p);
            if (!x(c.data(k)) && !x(c.data(k).secondary)) {
                o++, c.css(c.data(k));
                var C = c.data(k).secondary;
                setTimeout(function() {
                    c.bind(p, h).css(C)
                })
            } else
                d.queue = !1;
            x(q) || (o++, b.apply(c, [q, {duration: d.duration, easing: jQuery.easing[d.easing] ? d.easing : jQuery.easing.swing ? "swing" : "linear", complete: r, queue: d.queue}]));
            return!0
        })
    }, jQuery.fn.animate.defaults = {}, jQuery.fn.stop = function(a, b, d) {
        if (!q)
            return c.apply(this, [a, b]);
        a && this.queue([]);
        var e = {};
        for (var g = f.length - 1; g >= 0; g--)
            e[f[g] + "transition-property"] = "none", e[f[g] + "transition-duration"] = "", e[f[g] + "transition-timing-function"] = "";
        this.each(function() {
            var g = jQuery(this), h = window.getComputedStyle(this, null), j = {}, l;
            if (!x(g.data(k)) && !x(g.data(k).secondary)) {
                var m = g.data(k);
                if (b) {
                    j = m.secondary;
                    if (!d && typeof m.meta.left_o !== undefined || typeof m.meta.top_o !== undefined) {
                        j.left = typeof m.meta.left_o !== undefined ? m.meta.left_o : "auto", j.top = typeof m.meta.top_o !== undefined ? m.meta.top_o : "auto";
                        for (l = f.length - 1; l >= 0; l--)
                            j[f[l] + "transform"] = ""
                    }
                } else
                    for (var n in g.data(k).secondary) {
                        n = n.replace(i, "-$1").toLowerCase(), j[n] = h.getPropertyValue(n);
                        if (!d && /matrix/i.test(j[n])) {
                            var o = j[n].replace(/^matrix\(/i, "").split(/, |\)$/g);
                            j.left = parseFloat(o[4]) + parseFloat(g.css("left")) + "px" || "auto", j.top = parseFloat(o[5]) + parseFloat(g.css("top")) + "px" || "auto";
                            for (l = f.length - 1; l >= 0; l--)
                                j[f[l] + "transform"] = ""
                        }
                    }
                g.unbind(p).css(e).css(j).data(k, null)
            } else
                c.apply(g, [a, b])
        });
        return this
    }
})(jQuery, jQuery.fn.animate, jQuery.fn.stop)
function supports_html5_storage() {
    try {
        return'localStorage'in window && window['localStorage'] !== null;
    } catch (e) {
        return false;
    }
}
function findChildNavNodes() {
    $('.nb .m li a[rel=next]').each(function() {
        $(this).after("<span class='f c'><img alt='' src='/img/sui.gif'/></span>");
    });
}
function slideNav(offset) {
    $('html,body').animate({scrollTop: 0}, 0);
    $('#navcon').animate({left: offset, useTranslate3d: true}, 333, function() {
        if (slideBack) {
            $(this).children(".s").find('div:last-child').remove();
        } else {
            $(this).children(".s").find('div:first-child').remove();
        }
        $('#navcon').css({"left": 0});
        $('#navcon').children().find(".nb").css({"left": 0});
    });
}
findChildNavNodes();
var nav = $(document.getElementById('nav'));
if (nav[0]) {
    var navUrl = $("#nav").attr("data-id")
    var navId = "nav-" + navUrl.replace(/[a-z-./\s]/g, "");
    var slideBack
    var clickEventType = Modernizr.touch ? 'touchend' : 'click';
    if ($('html').hasClass('ios5 smallscreen')) {
        clickEventType = 'click';
    }
    if ($('html').hasClass('smallscreen')) {
        nav.appendTo("body");
    }
    nav.on('' + clickEventType + '', '.h li:last-child a', function(e) {
        e.preventDefault();
        window.location = $(this).attr("href");
    });
    nav.on(clickEventType, '.h li:not(:last-child)', function(e) {
        e.preventDefault();
        slideBack = true
        $(this).find("a").removeAttr("href");
        var currentNavBlock = $(this).parents("div");
        var currentNavBlockId = $(this).parents("div").attr("id");
        var selectedId = $(this).find("a").attr("id").replace("b", "l");
        var offset = currentNavBlock.position().left - currentNavBlock.width();
        if ($('#n' + selectedId).length <= 0) {
            if (supports_html5_storage() && localStorage.getItem("current-nav") == navId) {
                $("#" + currentNavBlockId).before('<div id="before" class="nb" style="position:absolute;left:' + offset + 'px" />');
                $("#before").html($(localStorage.getItem(navId)).find('#n' + selectedId).html()).attr("id", "n" + selectedId);
                nav.find('.s').css({"height": $('#n' + selectedId).height()});
                findChildNavNodes();
                if (currentNavBlock.position().left == 0) {
                    slideNav(Math.abs(offset));
                } else {
                    slideNav("-" + offset);
                }
            } else {
                $.ajax({url: navUrl, beforeSend: function() {
                        $("#" + selectedId.replace("l", "b")).parent("li").addClass("spin");
                        return;
                    }, success: function(html) {
                        $("#" + currentNavBlockId).before('<div id="before" class="nb" style="position:absolute;left:' + offset + 'px" />');
                        $("#before").html($(html).find('#n' + selectedId).html()).attr("id", "n" + selectedId);
                        nav.find('.s').css({"height": $('#n' + selectedId).height()});
                        findChildNavNodes();
                        if (currentNavBlock.position().left == 0) {
                            slideNav(Math.abs(offset));
                        } else {
                            slideNav("-" + offset);
                        }
                        if (supports_html5_storage()) {
                            localStorage.removeItem(localStorage.getItem("current-nav"));
                            localStorage.setItem(navId, html.replace(/\s+/g, " "));
                            localStorage.setItem("current-nav", navId);
                        }
                    }});
            }
            ;
        }
    });
    nav.on(clickEventType, '.nb ul li .c', function(e) {
        slideBack = false
        e.preventDefault();
        var currentNavBlock = $(this).parents("div");
        var currentNavBlockId = $(this).parents("div").attr("id");
        var selectedId = $(this).parent().find("a").attr("Id");
        var offset = currentNavBlock.position().left + currentNavBlock.width();
        if ($('#n' + selectedId).length <= 0) {
            if (supports_html5_storage() && localStorage.getItem("current-nav") == navId) {
                $("#" + currentNavBlockId).after('<div id="next" class="nb" style="position:absolute;left:' + offset + 'px" />');
                $("#next").html($(localStorage.getItem(navId)).find('#n' + selectedId).html()).attr("id", "n" + selectedId);
                findChildNavNodes();
                slideNav("-" + offset);
                nav.find('.s').css({"height": $('#n' + selectedId).height()});
            } else {
                $.ajax({url: navUrl, beforeSend: function() {
                        $("#" + selectedId).parents("li").addClass("spin");
                    }, success: function(html) {
                        $("#" + currentNavBlockId).after('<div id="next" class="nb" style="position:absolute; left:' + offset + 'px" />');
                        $("#next").html($(html).find('#n' + selectedId).html()).attr("id", "n" + selectedId);
                        nav.find('.s').css({"height": $('#n' + selectedId).height()});
                        findChildNavNodes();
                        slideNav("-" + offset);
                        if (supports_html5_storage()) {
                            localStorage.removeItem(localStorage.getItem("current-nav"));
                            localStorage.setItem(navId, html.replace(/\s+/g, " "));
                            localStorage.setItem("current-nav", navId);
                        }
                    }});
            }
            ;
        }
    });
    var mobsearch = $('#search');
    var mobsearchx = mobsearch.find('.ix');
    var con = $("#container");
    var wrap = $("#page-wrap");
    $(document).on(clickEventType, '#nav-mobile-open', function(e) {
        if (con.hasClass("on")) {
            wrap.removeClass("hdn");
            nav.removeClass("show");
            con.removeClass("on");
        } else {
            wrap.addClass("hdn");
            nav.addClass("show");
            con.addClass("on");
            mobsearchx.hide();
            mobsearch.removeClass('on');
        }
        e.preventDefault();
    });
}
var JSON;
JSON || (JSON = {});
(function() {
    function k(a) {
        return a < 10 ? "0" + a : a
    }
    function o(a) {
        p.lastIndex = 0;
        return p.test(a) ? '"' + a.replace(p, function(a) {
            var c = r[a];
            return typeof c === "string" ? c : "\\u" + ("0000" + a.charCodeAt(0).toString(16)).slice(-4)
        }) + '"' : '"' + a + '"'
    }
    function l(a, j) {
        var c, d, h, m, g = e, f, b = j[a];
        b && typeof b === "object" && typeof b.toJSON === "function" && (b = b.toJSON(a));
        typeof i === "function" && (b = i.call(j, a, b));
        switch (typeof b) {
            case"string":
                return o(b);
            case"number":
                return isFinite(b) ? String(b) : "null";
            case"boolean":
            case"null":
                return String(b);
            case"object":
                if (!b)
                    return"null";
                e += n;
                f = [];
                if (Object.prototype.toString.apply(b) === "[object Array]") {
                    m = b.length;
                    for (c = 0; c < m; c += 1)
                        f[c] = l(c, b) || "null";
                    h = f.length === 0 ? "[]" : e ? "[\n" + e + f.join(",\n" + e) + "\n" + g + "]" : "[" + f.join(",") + "]";
                    e = g;
                    return h
                }
                if (i && typeof i === "object") {
                    m = i.length;
                    for (c = 0; c < m; c += 1)
                        typeof i[c] === "string" && (d = i[c], (h = l(d, b)) && f.push(o(d) + (e ? ": " : ":") + h))
                } else
                    for (d in b)
                        Object.prototype.hasOwnProperty.call(b, d) && (h = l(d, b)) && f.push(o(d) + (e ? ": " : ":") + h);
                h = f.length === 0 ? "{}" : e ? "{\n" + e + f.join(",\n" + e) + "\n" + g + "}" : "{" + f.join(",") + "}";
                e = g;
                return h
            }
    }
    if (typeof Date.prototype.toJSON !== "function")
        Date.prototype.toJSON = function() {
            return isFinite(this.valueOf()) ? this.getUTCFullYear() + "-" + k(this.getUTCMonth() + 1) + "-" + k(this.getUTCDate()) + "T" + k(this.getUTCHours()) + ":" + k(this.getUTCMinutes()) + ":" + k(this.getUTCSeconds()) + "Z" : null
        }, String.prototype.toJSON = Number.prototype.toJSON = Boolean.prototype.toJSON = function() {
            return this.valueOf()
        };
    var q = /[\u0000\u00ad\u0600-\u0604\u070f\u17b4\u17b5\u200c-\u200f\u2028-\u202f\u2060-\u206f\ufeff\ufff0-\uffff]/g, p = /[\\\"\x00-\x1f\x7f-\x9f\u00ad\u0600-\u0604\u070f\u17b4\u17b5\u200c-\u200f\u2028-\u202f\u2060-\u206f\ufeff\ufff0-\uffff]/g, e, n, r = {"\u0008": "\\b", "\t": "\\t", "\n": "\\n", "\u000c": "\\f", "\r": "\\r", '"': '\\"', "\\": "\\\\"}, i;
    if (typeof JSON.stringify !== "function")
        JSON.stringify = function(a, j, c) {
            var d;
            n = e = "";
            if (typeof c === "number")
                for (d = 0; d < c; d += 1)
                    n += " ";
            else
                typeof c === "string" && (n = c);
            if ((i = j) && typeof j !== "function" && (typeof j !== "object" || typeof j.length !== "number"))
                throw Error("JSON.stringify");
            return l("", {"": a})
        };
    if (typeof JSON.parse !== "function")
        JSON.parse = function(a, e) {
            function c(a, d) {
                var g, f, b = a[d];
                if (b && typeof b === "object")
                    for (g in b)
                        Object.prototype.hasOwnProperty.call(b, g) && (f = c(b, g), f !== void 0 ? b[g] = f : delete b[g]);
                return e.call(a, d, b)
            }
            var d, a = String(a);
            q.lastIndex = 0;
            q.test(a) && (a = a.replace(q, function(a) {
                return"\\u" + ("0000" + a.charCodeAt(0).toString(16)).slice(-4)
            }));
            if (/^[\],:{}\s]*$/.test(a.replace(/\\(?:["\\\/bfnrt]|u[0-9a-fA-F]{4})/g, "@").replace(/"[^"\\\n\r]*"|true|false|null|-?\d+(?:\.\d*)?(?:[eE][+\-]?\d+)?/g, "]").replace(/(?:^|:|,)(?:\s*\[)+/g, "")))
                return d = eval("(" + a + ")"), typeof e === "function" ? c({"": d}, "") : d;
            throw new SyntaxError("JSON.parse");
        }
})();
(function() {
    var b = $$id("container") || document.body;
    $(b).delegate("a", "click", function() {
        if (this.rel.match('pop')) {
            var b = Math.floor(Math.random() * 1111111), a = this.rel.split("|"), c, d, e, f, g, h, i, j;
            switch (a[0].toLowerCase()) {
                case"live":
                    e = 1018;
                    d = 680;
                    g = (screen.height - 680) / 2;
                    h = (screen.width - 1018) / 2;
                    f = "0";
                    i = "yes";
                    j = "fullWindow";
                    c = screen.height < 770 ? "yes" : "no";
                    break;
                case"popup":
                    e = a[1], d = a[2], f = a[3] || "1", c = "no", i = "yes", g = (screen.height - a[2]) / 2, h = (screen.width - a[1]) / 2
            }
            window.open(this.href, b, "width=" + e + ",height=" + d + ",top=" + g + ",left=" + h + ",type=" + j + ",resizable=" + i + ",scrollbars=" + f + ",channelmode=" + c);
            return!1
        }
    })
})();
(function($) {
    var click = Modernizr.touch ? 'touchstart' : 'click';
    $.fn.bskyb_accordian = function(e) {
        return this.each(function() {
            var f = {acc: this, accx: $$c("acc-x", this), acchdr: $$c("acc-hdr", this), speed: 200, remember: !1, "static": !1}, oo = $.extend(f, $(this).data(), e), loc = window.localStorage;
            oo.acchdr[0].className += ' on';
            $(oo.acc).delegate('.acc-hdr', click, function(e) {
                e.stopImmediatePropagation();
                for (var b = oo.accx.length; --b >= 0; ) {
                    oo.accx[b].style.cssText = "height:0";
                }
                if (this.className.match(/\bon\b/)) {
                    this.className = this.className.replace(/\bon\b/g, '');
                }
                else {
                    for (var i = oo.acchdr.length; --i >= 0; ) {
                        oo.acchdr[i].className = oo.acchdr[i].className.replace(/\bon\b/g, '');
                    }
                    this.className += " on";
                    var b = $(this).next(".acc-x")[0], c = b.cloneNode(!0);
                    c.style.cssText = "height:auto;width:100%;position:absolute;visibility:hidden";
                    $$id("container").appendChild(c);
                    var e = c.offsetHeight;
                    c.parentNode.removeChild(c);
                    b.style.cssText = "height:auto;max-height:" + (e * 1.2) + "px"
                }
                (function() {
                    if (loc && oo.remember) {
                        for (var b = oo.acchdr.length, c = -1; --b >= 0; ) {
                            oo.acchdr[b].className.match(/\bon\b/g) && (c = b)
                        }
                        if (c <= 0) {
                            if (loc.removeItem(oo.remember), loc.jsstyles) {
                                loc.jsstyles = loc.jsstyles.replace("*" + oo.remember, ""), loc.jsstyles.match == "undefined" && loc.removeItem("jsstyles");
                            }
                        }
                        else {
                            loc[oo.remember] = c + "**<style id='" + oo.remember + "'>.acc[data-remember='" + oo.remember + "'] .acc-x {height:0px!important;}.acc[data-remember='" + oo.remember + "'] .acc-x:nth-child(" + (c * 2 + 2) + ") {height:auto!important}</style>", loc.jsstyles && !loc.jsstyles.match(oo.remember) ? loc.jsstyles += "*" + oo.remember : loc.jsstyles = "*" + oo.remember
                        }
                    }
                })()
            });
            (function() {
                if (loc && loc[oo.remember]) {
                    for (var b = loc[oo.remember].split("**")[0], c = oo.accx.length; --c >= 0; ) {
                        oo.acchdr[c].className = oo.acchdr[c].className.replace(/\bon\b/g, ""), oo.accx[c].style.cssText = "height:0px"
                    }
                    b != -1 && (oo.acchdr[b].className += " on", oo.accx[b].style.cssText = "height:auto");
                    (b = $$id(oo.remember)) && b.parentNode.removeChild(b)
                }
            })()
        })
    }
})(jQuery);
$($$c("enable-accordian", document)).bskyb_accordian();
(function($) {
    var loc = window.localStorage;
    function tabclick(click, x, oo, $this) {
        x != "x" && clearInterval(oo.intv);
        var a = $$t('a', $this)[0], ahref = a.href.split('#')[1];
        if (ahref) {
            var toshow = document.getElementById(ahref);
            click.preventDefault();
            $($$t('li', oo.anchors)).removeClass("on");
            $this.className += " on";
            $(oo.sections).addClass('hdn');
            $(toshow).removeClass('hdn');
            if ($this.getAttribute("data-load") && !$this.getAttribute("data-done")) {
                var aj = document.createElement('div');
                aj.className = 'ajax-loader';
                toshow.appendChild(aj);
                $.when($.ajax($this.getAttribute("data-load"))).done(function(data) {
                    $this.setAttribute("data-done", "done");
                    toshow.innerHTML = data;
                    aj.parentNode && aj.parentNode.removeChild(aj)
                })
            }
            ;
            (function() {
                if (Modernizr.localstorage && oo.remember) {
                    for (var b = oo.lis.length, c = -1; --b >= 0; ) {
                        oo.lis[b].className.match(/\bon\b/g) && (c = b);
                    }
                    if (c <= 0) {
                        if (loc.removeItem(oo.remember), loc.jsstyles) {
                            loc.jsstyles = loc.jsstyles.replace("*" + oo.remember, ""), loc.jsstyles.match == "undefined" && loc.removeItem("jsstyles");
                        }
                    } else {
                        loc[oo.remember] = c + "**<style id='" + oo.remember + "'>.enable-tabs[data-remember='" + oo.remember + "'] .tab-x {display:none}.enable-tabs[data-remember='" + oo.remember + "'] .tab-x:nth-child(" + (c + 2) + ") {display:block}</style>", loc.jsstyles && !loc.jsstyles.match(oo.remember) ? loc.jsstyles += "*" + oo.remember : loc.jsstyles = "*" + oo.remember
                    }
                }
            })()
        }
    }
    $.fn.bskyb_tabs = function(options) {
        return this.each(function() {
            var defaults = {tabs: $(this), anchors: $$c('tab-ul', this)[0], sections: $$c('tab-x', this), remember: false, auto: false};
            var oo = $.extend(defaults, $(this).data(), options);
            oo.lis = $$t('li', oo.anchors);
            !$$c('on', oo.anchors)[0] && ($$t('li', oo.anchors)[0].className += " on");
            var clicktype = "click";
            $(oo.anchors).delegate("li", clicktype, function(click, x) {
                !this.className.match(/disabled/) && tabclick(click, x, oo, this)
            });
            (function() {
                if (Modernizr.localstorage && loc[oo.remember]) {
                    for (var b = loc[oo.remember].split("**")[0], c = oo.lis.length; --c >= 0; ) {
                        oo.lis[c].className = oo.lis[c].className.replace(/\bon\b/g, ""), oo.sections[c].style.cssText = "display:none"
                    }
                    b != -1 && (oo.lis[b].className += " on", oo.sections[b].style.cssText = "display:block");
                    $(oo.lis[b]).trigger("click");
                    (b = $$id(oo.remember)) && b.parentNode.removeChild(b)
                }
            })();
            if (oo.auto) {
            }
            (function() {
                if (location.hash) {
                    var as = $$c('a', oo.anchors), l = location.href.replace('#');
                    for (var i = as.length; --i >= 0; ) {
                        l.match(as[i].href.split('#')[1]) && $(as[i].parentNode).trigger('click');
                    }
                }
            })();
        })
    }
})(jQuery);
$($$c("enable-tabs", document)).bskyb_tabs();
(function() {
    window.dropdownmenuisopen = false;
    function findPos(a) {
        var b = 0, curtop = 0;
        if (a.offsetParent) {
            do
                b += a.offsetLeft, curtop += a.offsetTop;
            while (a = a.offsetParent);
            return[b, curtop]
        }
    }
    function sweep(sel) {
        for (var len = -1; ++len < sel.length; ) {
            sel[len].className && sel[len].className.match(/selected/) && (pos = len);
            sel[len].className = sel[len].className.replace(/selected/, '');
        }
    }
    function close(e) {
        window.dropdownmenuisopen = false;
        if (!e.target.className.match(/select-head|sui-ico-date/)) {
            var p = e.target;
            var i = 6;
            while (--i >= 0) {
                if (p && p.className && p.className.match(/\bselect\b|video-datesearch/)) {
                    break;
                } else {
                    if (p && i > 0) {
                        p = p.parentNode
                    } else {
                        for (i in selects) {
                            selects[i]['option'].style.cssText = "visibility:hidden;opacity:0";
                        }
                        $('.video-datesearch .ix,.video-datesearch .sui-ico-date').removeClass('on');
                        break;
                    }
                    ;
                }
            }
        }
    }
    function clickme(e) {
        window.dropdownmenuisopen = true;
        var disabled = $(e.target).closest('.select')[0].className.match(/disabled/);
        if (!disabled) {
            function hide() {
                sel['option'].style.cssText = "visisibity:hidden;opacity:0;filter:alpha(opacity:0)";
                $(window).unbind('scroll', hide);
                $(document.body).unbind('keydown');
            }
            $(document.body).unbind('keydown');
            for (i in selects) {
                selects[i]['option'].style.cssText = "visibility:hidden;opacity:0";
            }
            var width = this.offsetWidth;
            var id = this.id;
            var offset = document.body.scrollTop || document.documentElement.scrollTop;
            var sel = selects[id];
            var h = document.documentElement.clientHeight;
            var lh = sel['li'].length * sel['li'][1].offsetHeight;
            var seltop = sel['on'] && sel['on'].offsetTop;
            var top = offset + 12;
            var stop = findPos(sel['select'])[1] + 12;
            var left = findPos(sel['select'])[0] - (90 - (width / 2));
            var y = "auto";
            if (lh < h) {
                top = stop - (lh / 2);
                h = (lh + 46);
                top > offset ? top = top : top = offset + 12;
                selltop = 0;
                y = "hidden";
                top < offset + 12 && (top = offset + 12);
            }
            sel['option'].style.cssText = "left:" + left + "px;top:" + top + "px;visibility:visible;opacity:0.98;filter:alpha(opacity:98);height:" + (h - 24) + "px";
            sel['ul'].style.cssText = "height:" + (h - 46) + "px;overflow-y:" + y;
            sel['ul'].offsetTop;
            sel['ul'].focus();
            sel['on'] && (sel['ul'].scrollTop = ((seltop) - (stop - offset)) + 32);
            sel['option'].style.overflow = 'hidden';
            $(document.body).keydown(function(e) {
                if (e.keyCode === 40) {
                    sweep(sel['li']);
                    e.preventDefault();
                    pos = pos === sel['li'].length - 1 ? sel['li'].length - 1 : pos + 1;
                    sel['li'][pos].className += ' selected';
                    sel['on'] = sel['li'][pos];
                    sel['ul'].scrollTop = (sel['li'][pos].offsetTop)
                }
                if (e.keyCode === 38) {
                    sweep(sel['li']);
                    e.preventDefault();
                    pos = pos === 1 ? 1 : pos - 1;
                    sel['li'][pos].className += ' selected';
                    sel['on'] = sel['li'][pos];
                    sel['ul'].scrollTop = (sel['li'][pos].offsetTop)
                }
                if (e.keyCode === 13) {
                    var ahref = sel['on'].getElementsByTagName('a');
                    ahref[0] && (location.href = ahref[0].href)
                }
                ;
                e.keyCode === 27 && hide();
            })
            $(window).scroll($.throttle(250, hide));
        }
    }
    var pos = 0, selects = {}, sel = $('.select'), b = document.body;
    for (var i = sel.length; --i >= 0; ) {
        if (!Modernizr.touch || sel[i].className.match(/\bt2\b/)) {
            var id = sel[i].id;
            selects[id] = {};
            var opts = $$t('option', sel[i]);
            var select = document.createElement('div');
            selects[id].select = select;
            selects[id].select.className = sel[i].className;
            sel[i].disabled && (selects[id].select.className += ' disabled');
            selects[id].select.id = id;
            var option = document.createElement('div');
            selects[id].option = option;
            option.className = 'option';
            option.id = sel[i].id + "-option";
            option.style.height = "5px";
            var ul = document.createElement('ul');
            selects[id].ul = ul;
            var html = "";
            var selected = false;
            for (var j = -1; ++j < opts.length; ) {
                var li = opts[j], og = '';
                if (opts[j].parentNode.nodeName.match(/optgroup/i) && opts[j] == opts[j].parentNode.getElementsByTagName('option')[0]) {
                    og = '<li class="optgroup">' + opts[j].parentNode.label + '</li>';
                }
                var c = opts[j].hasAttribute && opts[j].hasAttribute('selected') ? "selected " + opts[j].className : opts[j].className;
                var d = opts[j].disabled ? "disabled " + opts[j].className : "";
                var h2h = opts[j].getAttribute('data-h2h');
                h2h = h2h ? "data-h2h=" + h2h : "";
                html += og + "<li " + h2h + " class='" + c + " " + d + "'><a href='" + opts[j].value + "'>" + opts[j].innerHTML + "</a><img src='/img/sui.gif' class='sui sui-select'/></li>"
            }
            if (sel[i].getAttribute('data-inpage')) {
                var hidden = document.createElement('input');
                hidden.type = 'hidden';
                hidden.name = sel[i].id;
                sel[i].parentNode.insertBefore(hidden, sel[i]);
                $(option).delegate('a', 'click', selects[sel[i].id], function(e) {
                    e.preventDefault();
                    if (!this.parentNode.className.match(/disabled/)) {
                        var el = $$t('span', e.data.select)[0];
                        el.innerText ? el.innerText = this.innerText : el.textContent = this.textContent;
                        document.getElementsByName(e.data.select.id)[0].value = this.getAttribute('href');
                        e.data['option'].style.cssText = "visibility:hidden;opacity:0";
                        for (var k = e.data['li'].length; --k >= 0; ) {
                            if (e.data['li'][k].className.match(/selected/)) {
                                e.data['li'][k].className = e.data['li'][k].className.replace(/selected/, '');
                                break;
                            }
                        }
                        !this.parentNode.className.match(/selected/) && (this.parentNode.className += ' selected');
                    }
                });
            }
            ul.innerHTML = html;
            option.appendChild(ul);
            b.appendChild(option);
            selects[id].li = $$t('li', ul);
            selects[id].on = $$c('on', ul);
            var selected = false;
            if (sel[i].getAttribute('data-showselected')) {
                selected = $$c('selected', option)[0];
                selected = selected.innerText ? selected.innerText : selected.textContent;
            }
            select.innerHTML = "<h4 class='select-head'><span>" + (selected || sel[i].title) + " </span><img src='/img/sui.gif' class='sui sui-select' /></h4>";
            sel[i].parentNode.insertBefore(select, sel[i]);
            sel[i].parentNode.removeChild(sel[i]);
        }
        $(document).delegate('.select', 'click', clickme);
        $(document.body).click(close)
    }
})();
;
(function($)
{
    window.$call = $(['call']);
    $.event.special.load = {add: function(hollaback)
        {
            if (this.nodeType === 1 && this.tagName.toLowerCase() === 'img' && this.src !== '')
            {
                if (this.complete || this.readyState === 4)
                {
                    hollaback.handler.apply(this);
                }
                else if (this.readyState === 'uninitialized' && this.src.indexOf('data:') === 0)
                {
                    $(this).trigger('error');
                }
                else
                {
                    $(this).bind('load', hollaback.handler);
                }
            }
        }};
    $.fn.bskyb_thickbox = function(options)
    {
        var init = function(oo)
        {
            !document.getElementById('v5-thick') && $(document.getElementsByTagName('body')[0]).prepend('<div id="v5-thick"></div><div id="v5-thick-con" class="' + oo.addclass + '"><div id="v5-thick-inner"><div id="v5-thick-frame"><div id="v5-thick-content"></div><div id="v5-thick-caption"></div></div></div></div>');
            document.getElementById('v5-thick').offsetTop;
            $('#v5-thick').css({opacity: 0.6})
            function addCaption()
            {
                var c = document.getElementById('v5-thick-caption');
                c.style.cssText = "height:30px";
                oo.caption && (c.innerHTML = "<p>" + oo.caption + "</p>");
            }
            function sizeFrame(oo, f)
            {
                $('#v5-thick-content').css({width: oo.width, 'height': oo.height})
                $('#v5-thick-frame').css({width: oo.width, 'min-height': oo.height})
                $('#v5-thick-frame').css({opacity: 1});
                try
                {
                    $('#v5-thick-frame')[0].addEventListener('webkitTransitionEnd', function(event)
                    {
                        $('#v5-thick-content').css({opacity: 1});
                        addCaption();
                    }, false);
                    $('#v5-thick-frame')[0].addEventListener('transitionend', function(event)
                    {
                        $('#v5-thick-content').css({opacity: 1});
                        addCaption();
                    }, false);
                }
                catch (e)
                {
                    $('#v5-thick-content').animate({opacity: 1})
                    addCaption();
                }
                f && f();
            }
            var type = 'ajax', url = oo.load.split('?')[0];
            url.match(/^\/|^http/) && (type = 'ajax');
            url.match(/\.jpg$|\.jpeg$|\.gif$|\.png$/) && (type = 'img');
            url.match(/^#/) && (type = 'inline');
            oo.iframe && (type = 'iframe');
            switch (type)
            {
                case'img':
                    var im = new Image();
                    $(im).bind('load', function()
                    {
                        if (!oo.width)
                        {
                            oo.width = this.width
                        }
                        ;
                        if (!oo.height)
                        {
                            oo.height = this.height
                        }
                        if (oo.width > $(window).width())
                        {
                            oo.width = $(window).width() - 60
                        }
                        if (oo.height > $(window).height())
                        {
                            oo.height = $(window).height() - 60
                        }
                        $('#v5-thick-content').html($(this))
                        sizeFrame(oo, function()
                        {
                        })
                    }).attr({src: oo.load})
                    break;
                case'ajax':
                    if (!oo.width)
                    {
                        oo.width = 320
                    }
                    ;
                    $('#v5-thick-content').load(oo.load, function()
                    {
                        sizeFrame(oo, function()
                        {
                        })
                    });
                    break;
                case'inline':
                    if (!oo.width && !oo.height)
                    {
                        oo.width = 320;
                        oo.height = 200
                    }
                    ;
                    $('#v5-thick-content').html($(oo.load).html());
                    sizeFrame(oo, function()
                    {
                    })
                    break;
                case'iframe':
                    if (!oo.width && !oo.height)
                    {
                        oo.width = 320;
                        oo.height = 200
                    }
                    ;
                    $('#v5-thick-content').html("<iframe width='" + oo.width + "' height='" + oo.height + "' src='" + oo.load + "' FRAMEBORDER='0' SCROLLING='0'></iframe>");
                    sizeFrame(oo, function()
                    {
                    })
                    break;
                case'alert':
                    $('#v5-thick-content').html("<div class='bskyb-thickbox-alertbody'>" + oo.caption + "</div><button>" + oo.button + "</button>");
                    sizeFrame(oo, function()
                    {
                        $('#v5-thick-frame').css({opacity: 1});
                        $('#v5-thick-frame')[0].addEventListener('webkitTransitionEnd', function(event)
                        {
                            $('#v5-thick-content').css({opacity: 1})
                        }, false);
                    })
                    break;
            }
            this.exit = function()
            {
                $('#v5-thick-con').remove();
                $('#v5-thick').animate({opacity: 0}, 500, function() {
                    $(this).remove();
                });
                oo.width = false, oo.height = false;
                try
                {
                    if (!!Modernizr.csstransitions)
                    {
                        $('#v5-thick')[0].addEventListener('webkitTransitionEnd', function(event)
                        {
                            $('#v5-thick').remove();
                        }, false);
                    }
                    else
                    {
                        throw('notransition')
                    }
                }
                catch (e)
                {
                    $('#v5-thick').remove();
                }
            }
            $.fn.bskyb_thickbox.exit = this.exit;
            $('#v5-thick-con').click(function(e)
            {
                e.preventDefault();
                e.target.id == 'v5-thick-inner' && exit();
            })
            $(document).keyup(function(e)
            {
                e == null ? keycode = event.keyCode : keycode = e.which
                keycode == 27 && exit();
            });
        }
        return this.each(function()
        {
            var defaults = {el: $(this), width: false, height: false, iframe: false, type: false, load: $(this).attr('href'), caption: $(this).attr('title'), scroll: false, button: 'close', showclose: false, addclass: '', callback: false};
            $.extend(defaults, $(this).data(), options);
            defaults.el[0] == 'call' ? init(defaults) : defaults.el.click(function(e)
            {
                e.preventDefault();
                init(defaults)
            });
        })
    }
})(jQuery)
$($$c("enable-thickbox", document)).bskyb_thickbox();
(function(c) {
    c.fn.extend({autocomplete: function(a, b) {
            var p = typeof a == "string", b = c.extend({}, c.Autocompleter.defaults, {url: p ? a : null, data: p ? null : a, delay: p ? c.Autocompleter.defaults.delay : 10, max: b && !b.scroll ? 10 : 150}, b);
            b.highlight = b.highlight || function(a) {
                return a
            };
            b.formatMatch = b.formatMatch || b.formatItem;
            return this.each(function() {
                new c.Autocompleter(this, b)
            })
        }, result: function(a) {
            return this.bind("result", a);
        }, search: function(a) {
            return this.trigger("search", [a])
        }, flushCache: function() {
            return this.trigger("flushCache")
        }, setOptions: function(a) {
            return this.trigger("setOptions", [a])
        }, unautocomplete: function() {
            return this.trigger("unautocomplete")
        }});
    c.Autocompleter = function(a, b) {
        function p() {
            var a = i.selected();
            if (!a)
                return false;
            var c = a.result;
            k = c;
            if (b.multiple) {
                var d = r(f.val());
                d.length > 1 && (c = d.slice(0, d.length - 1).join(b.multipleSeparator) + b.multipleSeparator + c);
                c += b.multipleSeparator
            }
            f.val(c);
            n();
            f.trigger("result", [a.data, a.value]);
            return true
        }
        function l(a, c) {
            if (s == d)
                i.hide();
            else {
                var g = f.val();
                if (c || g != k)
                    k = g, g = j(g), g.length >= b.minChars ? (f.addClass(b.loadingClass), b.matchCase || (g = g.toLowerCase()), e(g, h, n)) : (f.removeClass(b.loadingClass), i.hide())
            }
        }
        function r(a) {
            if (!a)
                return[""];
            var a = a.split(b.multipleSeparator), g = [];
            c.each(a, function(a, b) {
                c.trim(b) && (g[a] = c.trim(b))
            });
            return g
        }
        function j(a) {
            if (!b.multiple)
                return a;
            a = r(a);
            return a[a.length - 1]
        }
        function n() {
            var g = i.visible();
            i.hide();
            clearTimeout(o);
            f.removeClass(b.loadingClass);
            b.mustMatch && f.search(function(a) {
                a || (b.multiple ? (a = r(f.val()).slice(0, -1), f.val(a.join(b.multipleSeparator) + (a.length ? b.multipleSeparator : ""))) : f.val(""))
            });
            g && c.Autocompleter.Selection(a, a.value.length, a.value.length)
        }
        function h(d, e) {
            if (e && e.length && g) {
                f.removeClass(b.loadingClass);
                i.display(e, d);
                var h = e[0].value;
                b.autoFill && j(f.val()).toLowerCase() == d.toLowerCase() && s != q && (f.val(f.val() + h.substring(j(k).length)), c.Autocompleter.Selection(a, k.length, k.length + h.length));
                i.show()
            } else
                n()
        }
        function e(g, d, e) {
            b.matchCase || (g = g.toLowerCase());
            var h = m.load(g);
            if (h && h.length)
                d(g, h);
            else if (typeof b.url == "string" && b.url.length > 0) {
                var f = {timestamp: +new Date};
                c.each(b.extraParams, function(a, b) {
                    f[a] = typeof b == "function" ? b() : b
                });
                c.ajax({mode: "abort", port: "autocomplete" + a.name, dataType: b.dataType, url: b.url, data: c.extend({'term': j(g)}, f), success: function(a) {
                        if (a && a.terms) {
                            a = a.terms;
                            var temp = "";
                            for (var i = -1; ++i < a.length; ) {
                                temp += a[i].t + "\n";
                            }
                            a = temp;
                        }
                        else {
                            $('.' + b.resultsClass).hide().find('ul').html("");
                            return false
                        }
                        var e;
                        if (!(e = b.parse && b.parse(a))) {
                            e = [];
                            for (var a = a.split("\n"), h = 0; h < a.length; h++) {
                                var f = c.trim(a[h]);
                                f && (f = f.split("|"), e[e.length] = {data: f, value: f[0], result: b.formatResult && b.formatResult(f, f[0]) || f[0]})
                            }
                        }
                        m.add(g, e);
                        d(g, e);
                    }})
            } else
                i.emptyList(), e(g)
        }
        var d, q;
        d = 46;
        q = 8;
        var f = c(a).attr("autocomplete", "off").addClass(b.inputClass), o, k = "", m = c.Autocompleter.Cache(b), g = 0, s, u = {mouseDownOnSelect: false}, i = c.Autocompleter.Select(b, a, p, u), t;
        c.browser.opera && c(a.form).bind("submit.autocomplete", function() {
            if (t)
                return t = false
        });
        f.bind((c.browser.opera ? "keypress" : "keydown") + ".autocomplete", function(a) {
            s = a.keyCode;
            switch (a.keyCode) {
                case 38:
                    a.preventDefault();
                    i.visible() ? i.prev() : l(0, true);
                    break;
                case 40:
                    a.preventDefault();
                    i.visible() ? i.next() : l(0, true);
                    break;
                case 33:
                    a.preventDefault();
                    i.visible() ? i.pageUp() : l(0, true);
                    break;
                case 34:
                    a.preventDefault();
                    i.visible() ? i.pageDown() : l(0, true);
                    break;
                case b.multiple && c.trim(b.multipleSeparator) == "," && 188:
                case 9:
                case 13:
                    if (p())
                        return a.preventDefault(), t = true, false;
                    break;
                case 27:
                    i.hide();
                    break;
                default:
                    clearTimeout(o), o = setTimeout(l, b.delay)
                }
        }).focus(function() {
            g++
        }).blur(function() {
            g = 0;
            u.mouseDownOnSelect || (clearTimeout(o), o = setTimeout(n, 200))
        }).on('click', function() {
            g++ > 1 && !i.visible() && l(0, true)
        }).bind("search", function() {
            function a(g, c) {
                var e;
                if (c && c.length)
                    for (var d = 0; d < c.length; d++)
                        if (c[d].result.toLowerCase() == g.toLowerCase()) {
                            e = c[d];
                            break
                        }
                typeof b == "function" ? b(e) : f.trigger("result", e && [e.data, e.value])
            }
            var b = arguments.length > 1 ? arguments[1] : null;
            c.each(r(f.val()), function(b, g) {
                e(g, a, a)
            })
        }).bind("flushCache", function() {
            m.flush()
        }).bind("setOptions", function(a, g) {
            c.extend(b, g);
            "data"in g && m.populate()
        }).bind("unautocomplete", function() {
            i.unbind();
            f.unbind();
            c(a.form).unbind(".autocomplete")
        })
    };
    c.Autocompleter.defaults = {inputClass: "ac_input", resultsClass: "ac_results", loadingClass: "ac_loading", append: false, minChars: 1, delay: 200, matchCase: false, matchSubset: true, matchContains: false, submitOnClick: false, cacheLength: 10, max: 100, mustMatch: false, extraParams: {}, selectFirst: true, formatItem: function(a) {
            return a[0]
        }, formatMatch: null, autoFill: false, width: 0, multiple: false, multipleSeparator: ", ", highlight: function(a, b) {
            if (!a) {
                a = ''
            }
            ;
            return a.replace(RegExp("(?![^&;]+;)(?!<[^<>]*)(" + b.replace(/([\^\$\(\)\[\]\{\}\*\.\+\?\|\\])/gi, "\\$1") + ")(?![^<>]*>)(?![^&;]+;)", "gi"), "<strong>$1</strong>")
        }, scroll: true, scrollHeight: 180};
    c.Autocompleter.Cache = function(a) {
        function b(b, c) {
            a.matchCase || (b = b.toLowerCase());
            var d = b.indexOf(c);
            a.matchContains == "word" && (d = b.toLowerCase().search("\\b" + c.toLowerCase()));
            return d == -1 ? false : d == 0 || a.matchContains
        }
        function p(b, c) {
            n > a.cacheLength && r();
            j[b] || n++;
            j[b] = c
        }
        function l() {
            if (!a.data)
                return false;
            var b = {}, e = 0;
            if (!a.url)
                a.cacheLength = 1;
            b[""] = [];
            for (var d = 0, j = a.data.length; d < j; d++) {
                var f = a.data[d], f = typeof f == "string" ? [f] : f, o = a.formatMatch(f, d + 1, a.data.length);
                if (o !== false) {
                    var k = o.charAt(0).toLowerCase();
                    b[k] || (b[k] = []);
                    f = {value: o, data: f, result: a.formatResult && a.formatResult(f) || o};
                    b[k].push(f);
                    e++ < a.max && b[""].push(f)
                }
            }
            c.each(b, function(b, g) {
                a.cacheLength++;
                p(b, g)
            })
        }
        function r() {
            j = {};
            n = 0
        }
        var j = {}, n = 0;
        setTimeout(l, 25);
        return{flush: r, add: p, populate: l, load: function(h) {
                if (!a.cacheLength || !n)
                    return null;
                if (!a.url && a.matchContains) {
                    var e = [], d;
                    for (d in j)
                        if (d.length > 0) {
                            var l = j[d];
                            c.each(l, function(a, c) {
                                b(c.value, h) && e.push(c)
                            })
                        }
                    return e
                } else if (j[h])
                    return j[h];
                else if (a.matchSubset)
                    for (d = h.length - 1; d >= a.minChars; d--)
                        if (l = j[h.substr(0, d)])
                            return e = [], c.each(l, function(a, c) {
                                b(c.value, h) && (e[e.length] = c)
                            }), e;
                return null
            }}
    };
    c.Autocompleter.Select = function(a, b, p, l) {
        var el = {}
        if (a.append) {
            el = $(b).parent(this)
        } else {
            el = document.body
        }
        function r() {
            var click = Modernizr.touch ? 'touchstart' : 'click';
            var autosubmit = a.submitonClick;
            o && (k = c("<div class='" + a.resultsClass + "'></div>").hide().css("position", "absolute").appendTo(el), m = c("<ul/>").appendTo(k).mouseover(function(a) {
                j(a).nodeName && j(a).nodeName.toUpperCase() == "LI" && (d = c("li", m).removeClass(h).index(j(a)), c(j(a)).addClass(h))
            }).on('click', function(a) {
                b.value = a.target.innerText ? a.target.innerText : a.target.textContent;
                autosubmit ? ($(b).parents('form').trigger('submit'), p()) : p();
            }).mousedown(function() {
                l.mouseDownOnSelect = true
            }).mouseup(function() {
                l.mouseDownOnSelect = false
            }), a.width > 0 && k.css("width", a.width), o = false)
        }
        function j(a) {
            for (a = a.target; a && a.tagName != "LI"; )
                a = a.parentNode;
            return!a ? [] : a
        }
        function n(b) {
            e.slice(d, d + 1).removeClass(h);
            d += b;
            d < 0 ? d = e.size() - 1 : d >= e.size() && (d = 0);
            b = e.slice(d, d + 1).addClass(h);
            if (a.scroll) {
                var c = 0;
                e.slice(0, d).each(function() {
                    c += this.offsetHeight
                });
                c + b[0].offsetHeight - m.scrollTop() > m[0].clientHeight ? m.scrollTop(c + b[0].offsetHeight - m.innerHeight()) : c < m.scrollTop() && m.scrollTop(c)
            }
        }
        var h;
        h = "ac_over";
        var e, d = -1, q, f = "", o = true, k, m;
        return{display: function(b, l) {
                r();
                q = b;
                f = l;
                m.empty();
                for (var k = a.max && a.max < q.length ? a.max : q.length, i = 0; i < k; i++)
                    if (q[i]) {
                        var j = a.formatItem(q[i].data, i + 1, k, q[i].value, f);
                        j !== false && (j = c("<li/>").html(a.highlight(j, f)).addClass(i % 2 == 0 ? "ac_even" : "ac_odd").appendTo(m)[0], c.data(j, "ac_data", q[i]))
                    }
                e = m.find("li");
                a.selectFirst && (e.slice(0, 1).addClass(h), d = 0);
                c.fn.bgiframe && m.bgiframe();
            }, next: function() {
                n(1)
            }, prev: function() {
                n(-1)
            }, pageUp: function() {
                d != 0 && d - 8 < 0 ? n(-d) : n(-8)
            }, pageDown: function() {
                d != e.size() - 1 && d + 8 > e.size() ? n(e.size() - 1 - d) : n(8)
            }, hide: function() {
                k && k.hide();
                e && e.removeClass(h);
                d = -1
            }, visible: function() {
                return k && k.is(":visible")
            }, current: function() {
                return this.visible() && (e.filter("." + h)[0] || a.selectFirst && e[0])
            }, show: function() {
                var d = c(b).offset();
                var padding = 12
                var t = d.top + b.offsetHeight
                var l = d.left
                if (a.append) {
                    t = "auto"
                    l = "auto"
                }
                k.css({width: typeof a.width == "string" || a.width > 0 ? a.width : c(b).width() + padding, top: t, left: l}).show();
                if (a.scroll && (m.scrollTop(0), m.css({maxHeight: a.scrollHeight, overflow: "auto"}), c.browser.msie && typeof document.body.style.maxHeight === "undefined")) {
                    var f = 0;
                    e.each(function() {
                        f += this.offsetHeight
                    });
                    d = f > a.scrollHeight;
                    m.css("height", d ? a.scrollHeight : f);
                    d || e.width(m.width() - parseInt(e.css("padding-left")) - parseInt(e.css("padding-right")))
                }
            }, selected: function() {
                var a = e && e.filter("." + h).removeClass(h);
                return a && a.length && c.data(a[0], "ac_data")
            }, emptyList: function() {
                m && m.empty()
            }, unbind: function() {
                k && k.remove()
            }}
    };
    c.Autocompleter.Selection = function(a, b, c) {
        if (a.createTextRange) {
            var l = a.createTextRange();
            l.collapse(true);
            l.moveStart("character", b);
            l.moveEnd("character", c);
            l.select()
        } else if (a.setSelectionRange) {
            a.setSelectionRange(b, c);
        }
        else if (a.selectionStart) {
            a.selectionStart = b, a.selectionEnd = c;
        }
        a.focus();
    }
})(jQuery);
;
(function(b) {
    b.fn.dragsort = function(k) {
        var d = b.extend({}, b.fn.dragsort.defaults, k), g = [], a = null, j = null;
        this.selector && b("head").append("<style type='text/css'>" + (this.selector.split(",").join(" " + d.dragSelector + ",") + " " + d.dragSelector) + " { cursor: move; }</style>");
        this.each(function(k, i) {
            b(i).is("table") && b(i).children().size() == 1 && b(i).children().is("tbody") && (i = b(i).children().get(0));
            var m = {draggedItem: null, placeHolderItem: null, pos: null, offset: null, offsetLimit: null, scroll: null, container: i, init: function() {
                    b(this.container).attr("data-listIdx", k).mousedown(this.grabItem).find(d.dragSelector).css({cursor: "move"});
                    b(this.container).children(d.itemSelector).each(function(a) {
                        b(this).attr("data-itemIdx", a)
                    })
                }, grabItem: function(e) {
                    if (!(e.which != 1 || b(e.target).is(d.dragSelectorExclude))) {
                        for (var c = e.target; !b(c).is("[data-listIdx='" + b(this).attr("data-listIdx") + "'] " + d.dragSelector); ) {
                            if (c == this)
                                return;
                            c = c.parentNode
                        }
                        a != null && a.draggedItem != null && a.dropItem();
                        b(e.target).css("cursor", "move");
                        a = g[b(this).attr("data-listIdx")];
                        a.draggedItem = b(c).closest(d.itemSelector);
                        var c = parseInt(a.draggedItem.css("marginTop")), f = parseInt(a.draggedItem.css("marginLeft"));
                        a.offset = a.draggedItem.offset();
                        a.offset.top = e.pageY - a.offset.top + (isNaN(c) ? 0 : c) - 1;
                        a.offset.left = e.pageX - a.offset.left + (isNaN(f) ? 0 : f) - 1;
                        if (!d.dragBetween)
                            c = b(a.container).outerHeight() == 0 ? Math.max(1, Math.round(0.5 + b(a.container).children(d.itemSelector).size() * a.draggedItem.outerWidth() / b(a.container).outerWidth())) * a.draggedItem.outerHeight() : b(a.container).outerHeight(), a.offsetLimit = b(a.container).offset(), a.offsetLimit.right = a.offsetLimit.left + b(a.container).outerWidth() - a.draggedItem.outerWidth(), a.offsetLimit.bottom = a.offsetLimit.top + c - a.draggedItem.outerHeight();
                        var c = a.draggedItem.height(), f = a.draggedItem.width(), h = a.draggedItem.attr("style");
                        a.draggedItem.attr("data-origStyle", h ? h : "");
                        d.itemSelector == "tr" ? (a.draggedItem.children().each(function() {
                            b(this).width(b(this).width())
                        }), a.placeHolderItem = a.draggedItem.clone().attr("data-placeHolder", true), a.draggedItem.after(a.placeHolderItem), a.placeHolderItem.children().each(function() {
                            b(this).css({borderWidth: 0, width: b(this).width() + 1, height: b(this).height() + 1}).html("&nbsp;")
                        })) : (a.draggedItem.after(d.placeHolderTemplate), a.placeHolderItem = a.draggedItem.next().css({height: c, width: f}).attr("data-placeHolder", true));
                        a.draggedItem.css({position: "absolute", opacity: 0.8, "z-index": 999, height: c, width: f});
                        b(g).each(function(a, b) {
                            b.createDropTargets();
                            b.buildPositionTable()
                        });
                        a.scroll = {moveX: 0, moveY: 0, maxX: b(document).width() - b(window).width(), maxY: b(document).height() - b(window).height()};
                        a.scroll.scrollY = window.setInterval(function() {
                            if (d.scrollContainer != window)
                                b(d.scrollContainer).scrollTop(b(d.scrollContainer).scrollTop() + a.scroll.moveY);
                            else {
                                var c = b(d.scrollContainer).scrollTop();
                                if (a.scroll.moveY > 0 && c < a.scroll.maxY || a.scroll.moveY < 0 && c > 0)
                                    b(d.scrollContainer).scrollTop(c + a.scroll.moveY), a.draggedItem.css("top", a.draggedItem.offset().top + a.scroll.moveY + 1)
                            }
                        }, 10);
                        a.scroll.scrollX = window.setInterval(function() {
                            if (d.scrollContainer != window)
                                b(d.scrollContainer).scrollLeft(b(d.scrollContainer).scrollLeft() + a.scroll.moveX);
                            else {
                                var c = b(d.scrollContainer).scrollLeft();
                                if (a.scroll.moveX > 0 && c < a.scroll.maxX || a.scroll.moveX < 0 && c > 0)
                                    b(d.scrollContainer).scrollLeft(c + a.scroll.moveX), a.draggedItem.css("left", a.draggedItem.offset().left + a.scroll.moveX + 1)
                            }
                        }, 10);
                        a.setPos(e.pageX, e.pageY);
                        b(document).bind("selectstart", a.stopBubble);
                        b(document).bind("mousemove", a.swapItems);
                        b(document).bind("mouseup", a.dropItem);
                        d.scrollContainer != window && b(window).bind("DOMMouseScroll mousewheel", a.wheel);
                        return false
                    }
                }, setPos: function(e, c) {
                    var f = c - this.offset.top, h = e - this.offset.left;
                    d.dragBetween || (f = Math.min(this.offsetLimit.bottom, Math.max(f, this.offsetLimit.top)), h = Math.min(this.offsetLimit.right, Math.max(h, this.offsetLimit.left)));
                    this.draggedItem.parents().each(function() {
                        if (b(this).css("position") != "static" && (!b.browser.mozilla || b(this).css("display") != "table")) {
                            var a = b(this).offset();
                            f -= a.top;
                            h -= a.left;
                            return false
                        }
                    });
                    if (d.scrollContainer == window)
                        c -= b(window).scrollTop(), e -= b(window).scrollLeft(), c = Math.max(0, c - b(window).height() + 5) + Math.min(0, c - 5), e = Math.max(0, e - b(window).width() + 5) +
                                Math.min(0, e - 5);
                    else
                        var l = b(d.scrollContainer), g = l.offset(), c = Math.max(0, c - l.height() - g.top) + Math.min(0, c - g.top), e = Math.max(0, e - l.width() - g.left) + Math.min(0, e - g.left);
                    a.scroll.moveX = e == 0 ? 0 : e * d.scrollSpeed / Math.abs(e);
                    a.scroll.moveY = c == 0 ? 0 : c * d.scrollSpeed / Math.abs(c);
                    this.draggedItem.css({top: f, left: h})
                }, wheel: function(e) {
                    if ((b.browser.safari || b.browser.mozilla) && a && d.scrollContainer != window) {
                        var c = b(d.scrollContainer), f = c.offset();
                        e.pageX > f.left && e.pageX < f.left + c.width() && e.pageY > f.top && e.pageY < f.top + c.height() && (f = e.detail ? e.detail * 5 : e.wheelDelta / -2, c.scrollTop(c.scrollTop() + f), e.preventDefault())
                    }
                }, buildPositionTable: function() {
                    var a = this.draggedItem == null ? null : this.draggedItem.get(0), c = [];
                    b(this.container).children(d.itemSelector).each(function(d, h) {
                        if (h != a) {
                            var g = b(h).offset();
                            g.right = g.left + b(h).width();
                            g.bottom = g.top + b(h).height();
                            g.elm = h;
                            c.push(g)
                        }
                    });
                    this.pos = c
                }, dropItem: function() {
                    if (a.draggedItem != null) {
                        b(a.container).find(d.dragSelector).css("cursor", "pointer");
                        a.placeHolderItem.before(a.draggedItem);
                        var e = a.draggedItem.attr("data-origStyle");
                        a.draggedItem.attr("style", e);
                        e == "" && a.draggedItem.removeAttr("style");
                        a.draggedItem.removeAttr("data-origStyle");
                        a.placeHolderItem.remove();
                        b("[data-dropTarget]").remove();
                        window.clearInterval(a.scroll.scrollY);
                        window.clearInterval(a.scroll.scrollX);
                        var c = false;
                        b(g).each(function() {
                            b(this.container).children(d.itemSelector).each(function(a) {
                                parseInt(b(this).attr("data-itemIdx")) != a && (c = true, b(this).attr("data-itemIdx", a))
                            })
                        });
                        c && d.dragEnd.apply(a.draggedItem);
                        a.draggedItem = null;
                        b(document).unbind("selectstart", a.stopBubble);
                        b(document).unbind("mousemove", a.swapItems);
                        b(document).unbind("mouseup", a.dropItem);
                        d.scrollContainer != window && b(window).unbind("DOMMouseScroll mousewheel", a.wheel);
                        return false
                    }
                }, stopBubble: function() {
                    return false
                }, swapItems: function(e) {
                    if (a.draggedItem == null)
                        return false;
                    a.setPos(e.pageX, e.pageY);
                    for (var c = a.findPos(e.pageX, e.pageY), f = a, h = 0; c == - 1 && d.dragBetween && h < g.length; h++)
                        c = g[h].findPos(e.pageX, e.pageY), f = g[h];
                    if (c == -1 || b(f.pos[c].elm).attr("data-placeHolder"))
                        return false;
                    j == null || j.top > a.draggedItem.offset().top || j.left > a.draggedItem.offset().left ? b(f.pos[c].elm).before(a.placeHolderItem) : b(f.pos[c].elm).after(a.placeHolderItem);
                    b(g).each(function(a, b) {
                        b.createDropTargets();
                        b.buildPositionTable()
                    });
                    j = a.draggedItem.offset();
                    return false
                }, findPos: function(a, b) {
                    for (var d = 0; d < this.pos.length; d++)
                        if (this.pos[d].left < a && this.pos[d].right > a && this.pos[d].top < b && this.pos[d].bottom > b)
                            return d;
                    return-1
                }, createDropTargets: function() {
                    d.dragBetween && b(g).each(function() {
                        var d = b(this.container).find("[data-placeHolder]"), c = b(this.container).find("[data-dropTarget]");
                        d.size() > 0 && c.size() > 0 ? c.remove() : d.size() == 0 && c.size() == 0 && (b(this.container).append(a.placeHolderItem.removeAttr("data-placeHolder").clone().attr("data-dropTarget", true)), a.placeHolderItem.attr("data-placeHolder", true))
                    })
                }};
            m.init();
            g.push(m)
        });
        return this
    };
    b.fn.dragsort.defaults = {itemSelector: "li", dragSelector: "li", dragSelectorExclude: "input, textarea, a[href]", dragEnd: function() {
        }, dragBetween: false, placeHolderTemplate: "<section>&nbsp;</section>", scrollContainer: window, scrollSpeed: 5}
})(jQuery);
(function($)
{
    var LEFT = "left", RIGHT = "right", UP = "up", DOWN = "down", NONE = "none", HORIZONTAL = "horizontal", VERTICAL = "vertical", AUTO = "auto", PHASE_START = "start", PHASE_MOVE = "move", PHASE_END = "end", PHASE_CANCEL = "cancel", SUPPORTS_TOUCH = 'ontouchstart'in window, START_EV = SUPPORTS_TOUCH ? 'touchstart' : 'mousedown', MOVE_EV = SUPPORTS_TOUCH ? 'touchmove' : 'mousemove', END_EV = SUPPORTS_TOUCH ? 'touchend' : 'mouseup', CANCEL_EV = 'touchcancel', PLUGIN_NS = 'TouchSwipe';
    var defaults = {fingers: 1, threshold: 75, maxTimeThreshold: null, swipe: null, swipeLeft: null, swipeRight: null, swipeUp: null, swipeDown: null, swipeStatus: null, click: null, triggerOnTouchEnd: true, allowPageScroll: "auto"};
    $.fn.swipe = function(method)
    {
        $this = $(this);
        var plugin = $this.data(PLUGIN_NS);
        if (plugin && typeof method === 'string')
        {
            if (plugin[method])
                return plugin[method].apply(this, Array.prototype.slice.call(arguments, 1));
            else
                $.error('Method ' + method + ' does not exist on jQuery.swipe');
        }
        else if (!plugin && (typeof method === 'object' || !method))
        {
            return init.apply(this, arguments);
        }
    }
    $.fn.swipe.defaults = defaults;
    $.fn.swipe.phases = {PHASE_START: PHASE_START, PHASE_MOVE: PHASE_MOVE, PHASE_END: PHASE_END, PHASE_CANCEL: PHASE_CANCEL}
    $.fn.swipe.directions = {LEFT: LEFT, RIGHT: RIGHT, UP: UP, DOWN: DOWN}
    $.fn.swipe.pageScroll = {NONE: NONE, HORIZONTAL: HORIZONTAL, VERTICAL: VERTICAL, AUTO: AUTO}
    function init(options)
    {
        if (options && (options.allowPageScroll == undefined && (options.swipe != undefined || options.swipeStatus != undefined)))
            options.allowPageScroll = NONE;
        if (!options)
            options = {};
        options = $.extend({}, $.fn.swipe.defaults, options);
        return this.each(function()
        {
            var $this = $(this);
            var plugin = $this.data(PLUGIN_NS);
            if (!plugin)
            {
                plugin = new TouchSwipe(this, options);
                $this.data(PLUGIN_NS, plugin);
            }
        });
    }
    function TouchSwipe(element, options)
    {
        var $element = $(element);
        var phase = "start";
        var triggerElementID = null;
        var fingerCount = 0;
        var start = {x: 0, y: 0};
        var end = {x: 0, y: 0};
        var delta = {x: 0, y: 0};
        var startTime = 0;
        var endTime = 0;
        var moveEvent;
        try
        {
            $element.bind(START_EV, touchStart);
            $element.bind(CANCEL_EV, touchCancel);
        }
        catch (e)
        {
            $.error('events not supported ' + START_EV + ',' + CANCEL_EV + ' on jQuery.swipe');
        }
        this.enable = function()
        {
            $element.bind(START_EV, touchStart);
            $element.bind(CANCEL_EV, touchCancel);
            return $element;
        }
        this.disable = function()
        {
            removeListeners();
            return $element;
        }
        this.destroy = function()
        {
            removeListeners();
            $element.data(PLUGIN_NS, null);
            return $element;
        }
        function touchStart(event)
        {
            event = event.originalEvent;
            var ret, evt = SUPPORTS_TOUCH ? event.touches[0] : event;
            phase = PHASE_START;
            if (SUPPORTS_TOUCH) {
                fingerCount = event.touches.length;
            }
            else
            {
                event.preventDefault();
            }
            distance = 0;
            direction = null;
            duration = 0;
            if (!SUPPORTS_TOUCH || fingerCount == options.fingers)
            {
                start.x = end.x = evt.pageX;
                start.y = end.y = evt.pageY;
                startTime = getTimeStamp();
                if (options.swipeStatus)
                    ret = triggerHandler(event, phase);
            }
            else
            {
                touchCancel(event);
            }
            if (ret === false)
            {
                phase = PHASE_CANCEL;
                triggerHandler(event, phase);
                return ret;
            }
            else
            {
                $element.bind(MOVE_EV, touchMove);
                $element.bind(END_EV, touchEnd);
            }
        }
        function touchMove(event)
        {
            event = event.originalEvent;
            if (phase == PHASE_END || phase == PHASE_CANCEL)
                return;
            var ret, evt = SUPPORTS_TOUCH ? event.touches[0] : event;
            end.x = evt.pageX;
            end.y = evt.pageY;
            endTime = getTimeStamp();
            direction = calculateDirection();
            if (SUPPORTS_TOUCH)
                fingerCount = event.touches.length;
            phase = PHASE_MOVE;
            validateDefaultEvent(event, direction);
            if (fingerCount == options.fingers || !SUPPORTS_TOUCH)
            {
                distance = calculateDistance();
                duration = calculateDuration();
                if (options.swipeStatus)
                    ret = triggerHandler(event, phase, direction, distance, duration);
                if (!options.triggerOnTouchEnd)
                {
                    var cancel = !validateSwipeTime();
                    if (validateSwipeDistance() === true)
                    {
                        phase = PHASE_END;
                        ret = triggerHandler(event, phase);
                    }
                    else if (cancel)
                    {
                        phase = PHASE_CANCEL;
                        triggerHandler(event, phase);
                    }
                }
            }
            else
            {
                phase = PHASE_CANCEL;
                triggerHandler(event, phase);
            }
            if (ret === false)
            {
                phase = PHASE_CANCEL;
                triggerHandler(event, phase);
            }
        }
        function touchEnd(event)
        {
            event = event.originalEvent;
            event.preventDefault();
            endTime = getTimeStamp();
            distance = calculateDistance();
            direction = calculateDirection();
            duration = calculateDuration();
            if (options.triggerOnTouchEnd || (options.triggerOnTouchEnd == false && phase == PHASE_MOVE))
            {
                phase = PHASE_END;
                if ((fingerCount == options.fingers || !SUPPORTS_TOUCH) && end.x != 0)
                {
                    var cancel = !validateSwipeTime();
                    if ((validateSwipeDistance() === true || validateSwipeDistance() === null) && !cancel)
                    {
                        triggerHandler(event, phase);
                    }
                    else if (cancel || validateSwipeDistance() === false)
                    {
                        phase = PHASE_CANCEL;
                        triggerHandler(event, phase);
                    }
                }
                else
                {
                    phase = PHASE_CANCEL;
                    triggerHandler(event, phase);
                }
            }
            else if (phase == PHASE_MOVE)
            {
                phase = PHASE_CANCEL;
                triggerHandler(event, phase);
            }
            $element.unbind(MOVE_EV, touchMove, false);
            $element.unbind(END_EV, touchEnd, false);
        }
        function touchCancel(event)
        {
            fingerCount = 0;
            start.x = 0;
            start.y = 0;
            end.x = 0;
            end.y = 0;
            delta.x = 0;
            delta.y = 0;
            endTime = 0;
            startTime = 0;
        }
        function triggerHandler(event, phase)
        {
            var ret;
            if (options.swipeStatus)
                ret = options.swipeStatus.call($element, event, phase, direction || null, distance || 0, duration || 0);
            if (phase == PHASE_CANCEL)
            {
                if (options.click && (fingerCount == 1 || !SUPPORTS_TOUCH) && (isNaN(distance) || distance == 0))
                    ret = options.click.call($element, event, event.target);
            }
            if (phase == PHASE_END)
            {
                if (options.swipe)
                {
                    ret = options.swipe.call($element, event, direction, distance, duration);
                }
                switch (direction)
                {
                    case LEFT:
                        if (options.swipeLeft)
                            ret = options.swipeLeft.call($element, event, direction, distance, duration);
                        break;
                    case RIGHT:
                        if (options.swipeRight)
                            ret = options.swipeRight.call($element, event, direction, distance, duration);
                        break;
                    case UP:
                        if (options.swipeUp)
                            ret = options.swipeUp.call($element, event, direction, distance, duration);
                        break;
                    case DOWN:
                        if (options.swipeDown)
                            ret = options.swipeDown.call($element, event, direction, distance, duration);
                        break;
                    }
            }
            if (phase == PHASE_CANCEL || phase == PHASE_END)
            {
                touchCancel(event);
            }
            if (ret !== undefined)
                return ret;
        }
        function validateSwipeDistance()
        {
            if (options.threshold !== null)
                return distance >= options.threshold;
            else
                return null;
        }
        function validateSwipeTime()
        {
            var result;
            if (options.maxTimeThreshold)
            {
                if (duration >= options.maxTimeThreshold)
                    result = false;
                else
                    result = true;
            }
            else
            {
                result = true;
            }
            return result;
        }
        function validateDefaultEvent(event, direction)
        {
            if (options.allowPageScroll == NONE)
            {
                event.preventDefault();
            }
            else
            {
                var auto = options.allowPageScroll == AUTO;
                switch (direction)
                {
                    case LEFT:
                        if ((options.swipeLeft && auto) || (!auto && options.allowPageScroll != HORIZONTAL))
                            event.preventDefault();
                        break;
                    case RIGHT:
                        if ((options.swipeRight && auto) || (!auto && options.allowPageScroll != HORIZONTAL))
                            event.preventDefault();
                        break;
                    case UP:
                        if ((options.swipeUp && auto) || (!auto && options.allowPageScroll != VERTICAL))
                            event.preventDefault();
                        break;
                    case DOWN:
                        if ((options.swipeDown && auto) || (!auto && options.allowPageScroll != VERTICAL))
                            event.preventDefault();
                        break;
                    }
            }
        }
        function calculateDuration()
        {
            return endTime - startTime;
        }
        function calculateDistance()
        {
            return Math.round(Math.sqrt(Math.pow(end.x - start.x, 2) + Math.pow(end.y - start.y, 2)));
        }
        function caluculateAngle()
        {
            var X = start.x - end.x;
            var Y = end.y - start.y;
            var r = Math.atan2(Y, X);
            var angle = Math.round(r * 180 / Math.PI);
            if (angle < 0)
                angle = 360 - Math.abs(angle);
            return angle;
        }
        function calculateDirection()
        {
            var angle = caluculateAngle();
            if ((angle <= 45) && (angle >= 0))
                return LEFT;
            else if ((angle <= 360) && (angle >= 315))
                return LEFT;
            else if ((angle >= 135) && (angle <= 225))
                return RIGHT;
            else if ((angle > 45) && (angle < 135))
                return DOWN;
            else
                return UP;
        }
        function getTimeStamp()
        {
            var now = new Date();
            return now.getTime();
        }
        function removeListeners()
        {
            $element.unbind(START_EV, touchStart);
            $element.unbind(CANCEL_EV, touchCancel);
            $element.unbind(MOVE_EV, touchMove);
            $element.unbind(END_EV, touchEnd);
        }}
}
)(jQuery);
var ls = {reorderTR: function(t) {
        function sortfunction(a, b) {
            var t1 = parseInt(a.getAttribute('data-pos').replace('T', ''), 10), t2 = parseInt(b.getAttribute('data-pos').replace('T', ''), 10)
            return t1 - t2;
        }
        var trs = $(t).find('tr');
        $.makeArray(trs);
        trs.sort(sortfunction);
        for (var i = -1; ++i < trs.length; ) {
            t[0].appendChild(trs[i])
        }
    }, init: function(p) {
        var click = Modernizr.touch ? 'touchstart' : 'click';
        ls.o = {};
        ls.o.trs = [];
        ls.o.fix = [];
        ls.o.game = {};
        ls.o.data = {};
        tables = function(table) {
            var t = $$t('tbody', $$id(table))[0];
            if (t) {
                var tr = $$t("tr", t);
                var i = tr.length;
                var h = $(tr[1]).outerHeight();
                for (var i = tr.length; --i >= 0; ) {
                    ls.o.trs.push(tr[i]);
                    if (tr[i].className.match(/ls-hr/)) {
                        var hr = document.createElement("hr");
                        hr.className = "ls-hr";
                        hr.style.cssText = "top:" + tr[i].offsetTop + "px";
                        tr[i].parentNode.parentNode.parentNode.appendChild(hr);
                    }
                    tr[i].setAttribute("h", h);
                    tr[i].setAttribute("rel", (i + 1));
                }
            }
        };
        betting = function(el) {
            $('#' + el).on('click', 'a', function(e) {
                var p = this.parentNode.parentNode;
                if (p && p.getAttribute('data-status') && p.getAttribute('data-status').match(/s/gi)) {
                    e.preventDefault();
                }
                if (this.getAttribute('data-status') && this.getAttribute('data-status').match(/s/gi)) {
                    e.preventDefault();
                }
            })
        };
        scores = function(el) {
            ls.o.list = $$id(el);
            ls.o.listsort = $$c('sort-con', ls.o.list);
            ls.o.userteams = $('#scores-custom').find('ul')[0];
            (function() {
                for (var i = -1; ++i < ls.o.listsort.length; ) {
                    var id = ls.o.listsort[i].id + '-order';
                    if (Modernizr.localstorage && window.localStorage[id]) {
                        var d = new Date, time = d.getMonth() + "-" + d.getDate(), order = JSON.parse(window.localStorage[id]);
                        if (order.date == time) {
                            order.order.reverse();
                            var style = $$id("ls-lge-order");
                            var frag = document.createDocumentFragment();
                            for (var j = order.order.length; --j >= 0; ) {
                                var d = $$id(order.order[j][0]);
                                d && frag.appendChild(d);
                            }
                            ls.o.listsort[i] && ls.o.listsort[i].appendChild(frag);
                            style && style.parentNode.removeChild(style);
                        } else {
                            window.localStorage.removeItem(id)
                        }
                    }
                }
            })()
            ls.o.fix = $$id("page") ? $$t('li', $$id("page")) : $$t('li', document);
            var l = window.localStorage;
            if (Modernizr.localstorage) {
                var d = new Date, time = {"date": "" + d.getMonth() + "-" + d.getDate() + "", "teams": []}, user = l.getItem("userscores") ? JSON.parse(l.getItem("userscores")) : time;
                l.setItem("userscores", JSON.stringify(user));
                user.date != time.date && (user = time, l.setItem("userscores", JSON.stringify(time)))
            }
            function clone() {
                if (user) {
                    var f = document.createDocumentFragment();
                    var i = user.teams.length;
                    while (--i >= 0) {
                        var b = $$id(user.teams[i]), c = b.cloneNode(true);
                        c.className += " minus";
                        c.id = c.id + "-clone";
                        f.appendChild(c);
                        $$t('div', $$t('a', b)[0])[3].className += " disabled";
                    }
                    ls.o.userteams && ls.o.userteams.insertBefore(f, ls.o.userteams.firstChild);
                    ls.o.userteams.style.cssText = "min-height:0"
                }
            }
            ls.o.userteams && clone();
            (function() {
                var state = false;
                try {
                    state = window.getComputedStyle(ls.o.list, null).getPropertyValue("min-width") == "1px"
                } catch (e) {
                }
                var b = $$id("bet-switch"), a = $$t('a', b);
                if (b) {
                    state ? (a[0].className = "fc", a[1].className = "on") : a[0].className = "fc on"
                }
                ;
                $(b).on('click', 'a', function(e) {
                    if (e.target.className.match(/fc/)) {
                        a[0].className = "fc on", a[1].className = "";
                        ls.o.list.className = "c1 off";
                        if (Modernizr.localstorage) {
                            window.localStorage["ls-betswitch"] = "off"
                        }
                    } else {
                        a[0].className = "fc", a[1].className = "on";
                        ls.o.list.className = "c1 on";
                        if (Modernizr.localstorage) {
                            window.localStorage["ls-betswitch"] = "on"
                        }
                    }
                    return false
                })
            })();
            $(ls.o.list).delegate(".ixt", 'click', function(e) {
                e.preventDefault();
                var li = this.parentNode;
                if (e.target.className.match('sui-odds') && e.target.nodeName.match(/a/i)) {
                    e.target.getAttribute('target') ? window.open(e.target.href) : location.href = e.target.href;
                    return false;
                }
                if (e.target.className.match(/\bixm\b/)) {
                    if (e.target.className.match(/disabled/) || li.className.match(/ft|\ba\b/)) {
                        return false
                    }
                    else {
                        if (li.className.match(/minus/g)) {
                            var t = $$t('div', $$t('a', $$id(li.id.replace("-clone", "")))[0])[3]
                            $(t).removeClass('disabled');
                            li.style.cssText = "opacity:0";
                            Modernizr.csstransitions ? li.parentNode.removeChild(li) : li.parentNode.removeChild(li);
                            user.teams.length == 1 && ($$c('note', ls.o.userteams)[0].style.display = "block");
                            if (l) {
                                var i = user.teams.length;
                                while (--i >= 0) {
                                    user.teams[i] == li.id.replace("-clone", "") && user.teams.splice(i, 1)
                                }
                                !user || l.setItem("userscores", JSON.stringify(user))
                            }
                        } else {
                            $$c('note', ls.o.userteams)[0].style.display = "none";
                            var c = li.cloneNode(true);
                            c.className += " minus";
                            c.style.cssText = "opacity:0";
                            c.id = li.id + "-clone";
                            e.target.className += " disabled";
                            ls.o.userteams.insertBefore(c, ls.o.userteams.firstChild);
                            c.offsetTop;
                            c.style.cssText = "opacity:1";
                            !user || user.teams.push(li.id), l.setItem("userscores", JSON.stringify(user))
                        }
                    }
                } else {
                    !li.className.match(/nolink/) && (location.href = li.getElementsByTagName("a")[0].href)
                }
            });
            $(ls.o.listsort).dragsort({itemSelector: ".ls-x:not(.ui-state-disabled)", dragSelector: '.ls-x:not(.ui-state-disabled) h4', dragEnd: function() {
                    var parent = $(this).parent(), id = parent.attr('id');
                    if (Modernizr.localstorage) {
                        var order = [];
                        var s = $$t('section', parent[0]), i = s.length;
                        while (--i >= 0) {
                            order[i] = [s[i].id, s[i].offsetTop]
                        }
                        var d = new Date, time = {"date": "" + d.getMonth() + "-" + d.getDate()};
                        time.order = order;
                        window.localStorage.setItem(id + '-order', JSON.stringify(time))
                    }
                }});
        };
        game = function(el) {
            ls.o.game.game = $$id("ls-match-hdr");
            ls.o.game.score = $$id("ls-match-score");
            ls.o.game.hs = $$id("ls-home-score");
            ls.o.game.as = $$id("ls-away-score");
            ls.o.game.hg = $$id("ls-home-goals");
            ls.o.game.ag = $$id("ls-away-goals");
            ls.o.game.att = $$id("ls-match-att");
            ls.o.game.status = $$id("ls-match-status");
            ls.o.game.ht = $$id("ls-match-ht");
            ls.o.game.desc = $$id("ls-match-desc");
            ls.o.game.players = $$c('ls-sq', document);
            $(".ls-match-squad").delegate(".ls-sq-exp,.sq-ix", "click", function(e) {
                if (e.target.className.match(/ls-sq-exp|sui/)) {
                    e.target.nodeName.match(/img/i) && (e.target = e.target.parentNode);
                    var form = $$c('ls-sq-form', this.parentNode.parentNode)[0];
                    if (form) {
                        if (!form.className.match(/form-/)) {
                            var f = e.target.innerText ? e.target.innerText.replace(/ /g, "") : e.target.textContent.replace(/ /g, "");
                            var ps = $$c('ls-sq', form.parentNode), df = document.createDocumentFragment(), len = form.className.match(/rugby/) ? 15 : 11;
                            for (var i = 0; i < len; i++) {
                                var d = document.createElement("div");
                                var css = ps[i].className.match(/\bon\b/) == "on" ? "on" : "";
                                d.className = css + " enable-tooltip sq-ix p-" + i;
                                d.id = "clone-" + ps[i].id;
                                d.setAttribute("href", $$t("a", ps[i])[0].href);
                                d.innerHTML = $$t('strong', ps[i])[0].innerHTML.replace(".", "");
                                d.setAttribute('data-text', ps[i].innerHTML);
                                df.appendChild(d);
                            }
                            var ievml = '<div class="l4"></div>';
                            if (bskyb.browser.match(/ie8|ie7/)) {
                                (document.namespaces.add('v', 'urn:schemas-microsoft-com:vml'));
                                ievml = '<v:oval style="position:absolute;top:84px;left:72px;width:58px;height:58px;z-index:2" strokeweight="2px" strokecolor="white" fillcolor="#6e9c48"/>'
                            }
                            ;
                            form.innerHTML = '<div class="pitch"><div class="l1"></div><div class="l2"> </div><div class="l3"></div><div class="l3 t2"></div>' + ievml + '<div class="l5"></div></div>';
                            form.appendChild(df);
                            form.offsetTop;
                            form.className += " form-" + f
                        }
                        e.target.className.match(/on/) ? (e.target.className = "ls-sq-exp", form.className = form.className.replace(/\b on\b/, "")) : (e.target.className += " on", form.className += " on")
                    }
                }
                if (e.target.className.match(/sq-ix/) || e.target.parentNode.className.match(/sq-ix/)) {
                    e.preventDefault();
                    var id = this.id || this.parentNode.id;
                    var parent = $(this).closest(".ls-match-squad");
                    id = id.replace("clone-", ""), el = $$id(id);
                    if (this.className.match(/\bon\b/) || $(this).closest("li").hasClass("on")) {
                        $$id(id + "-prof").className = "ls-sq-prof";
                        $("#clone-" + id).removeClass("on");
                        $("#" + id).removeClass("on")
                    }
                    else {
                        parent.find(".sq-ix,.ls-sq,.ls-sq-prof").removeClass("on");
                        $("#clone-" + id).addClass("on");
                        $("#" + id).addClass("on");
                        var exists = $$id(id + "-prof");
                        if (!exists) {
                            var d = document.createElement("div");
                            d.className = "ls-sq-prof ajax";
                            d.id = id + "-prof";
                            el.parentNode.insertBefore(d, el.nextSibling);
                            d.offsetTop;
                            d.className += " on";
                            $(d).load(this.getAttribute("href"), function() {
                                $(this).removeClass("ajax")
                            })
                        }
                        else {
                            exists.className.match(/on/) ? (exists.className = "ls-sq-prof", this.className = this.className.replace(/on/, "")) : exists.className += " on"
                        }
                    }
                }
            })
        };
        golfscores = function(el) {
            ls.o.golf = {};
            ls.o.golf.timestamp = $$id('time-stamp');
            ls.o.golf.tbody1 = $("#myscores-leaderboard").find("tbody"), ls.o.golf.tbody2 = $("#live-leaderboard").find("tbody");
            var trs = $("#live-leaderboard").find("tbody").find("tr"), trsX = $("tr");
            for (var i = -1; ++i < trs.length; ) {
                trs[i].setAttribute("rel", i * (bskyb.smallscreen ? 32 : 40));
                trs[i].setAttribute('data-pos', trs[i].getElementsByTagName('td')[0].innerHTML)
            }
            function addToLS(t) {
                var trs = t[0].getElementsByTagName("tr"), mygolfers = [];
                for (var i = trs.length; --i >= 0; ) {
                    mygolfers.push(trs[i].getAttribute("id").replace("-clone", ""))
                }
                Modernizr.localstorage && (window.localStorage.myleaderboard = mygolfers)
            }
            function getFromLS() {
                if (Modernizr.localstorage && window.localStorage.myleaderboard) {
                    var ids = window.localStorage.myleaderboard.split(","), mytable = $("#myscores-leaderboard"), tbody = mytable.find("tbody");
                    mytable.removeClass("hdn");
                    for (var i = ids.length; --i >= 0; ) {
                        var tr = $("#" + ids[i]);
                        var clone = tr.clone(true);
                        tr.addClass("on");
                        clone[0].id = clone[0].id + "-clone";
                        clone[0].style.cssText = "";
                        clone[0].className += " clone";
                        tbody.append(clone[0])
                    }
                    ls.reorderTR(tbody)
                }
            }
            getFromLS();
            $(".live-ldrb").on("click", ".sui-myfave", function() {
                var tr = this.parentNode.parentNode, mytable = $("#myscores-leaderboard"), tbody = mytable.find("tbody");
                if (!tr.className.match(/\bon\b/gi) && !tr.className.match(/\bclone\b/gi)) {
                    mytable.removeClass("hdn");
                    var clone = tr.cloneNode(true);
                    clone.id = clone.id + "-clone";
                    clone.style.cssText = "";
                    tbody.append(clone);
                    ls.reorderTR(tbody);
                    clone.className += " clone";
                    tr.className += " on"
                }
                if (tr.className.match(/\bclone\b/gi)) {
                    var id = tr.getAttribute("id").replace("-clone", "");
                    $("#" + id).removeClass("on");
                    tr.parentNode.removeChild(tr);
                    if (tbody.find("tr").length < 1) {
                        mytable.addClass("hdn")
                    }
                }
                if (tr.className.match(/\bon\b/gi)) {
                    console.log($("#" + tr.getAttribute("id") + "-clone"))
                }
                addToLS(tbody)
            });
        }
        tennisscores = function(el) {
            function cssblock(p, c) {
                if (p.styleSheet) {
                    p.styleSheet.cssText = c;
                }
                else {
                    p.innerHTML = c;
                }
            }
            ls.o.tennis = {};
            ls.o.tennis.con = $$id(el), ls.o.tennis.matches = $$c("sui-draw-fix", ls.o.tennis.con), ls.o.tennis.groups = $$c("sui-draw-fix-gp", ls.o.tennis.con), ls.o.tennis.style = document.createElement("style"), ls.o.tennis.style2 = document.createElement("style"), ls.o.tennis.filters = $(ls.o.tennis.con).find(".filters-radio");
            ls.o.tennis.style.id = "sui-draw-fix-styles", ls.o.tennis.style2.id = "sui-draw-fix-styles-hdr";
            ls.o.tennis.style.type = "text/css", ls.o.tennis.style2.type = "text/css";
            document.head.appendChild(ls.o.tennis.style);
            document.head.appendChild(ls.o.tennis.style2);
            $(ls.o.tennis.filters).on("click", "input", function(e) {
                var t = this.getAttribute("value");
                switch (t) {
                    case"tennisfaves":
                        cssblock(ls.o.tennis.style, "table[data-status='To Finish'],table[data-status='Finished'],table[data-status=NSY],table[data-status=inprogress]{display:none}");
                        Modernizr.localstorage && (window.sessionStorage.tennislivetabs = "3");
                        break;
                    case"Finished":
                        cssblock(ls.o.tennis.style, "table.tennisfave,table[data-status='To Finish'],table[data-status=NSY],table[data-status=inprogress]{display:none}");
                        Modernizr.localstorage && (window.sessionStorage.tennislivetabs = "2");
                        break;
                    case"inprogress":
                        cssblock(ls.o.tennis.style, "table[data-status=NSY],table[data-status=Finished]{display:none}");
                        Modernizr.localstorage && (window.sessionStorage.tennislivetabs = "1");
                        break;
                    case"all":
                        cssblock(ls.o.tennis.style, "table.tennisfave {display:none}");
                        Modernizr.localstorage && (window.sessionStorage.tennislivetabs = "0");
                        break;
                    case"NSY":
                        cssblock(ls.o.tennis.style, "table.tennisfave,table[data-status='To Finish'],table.tennisfave,table[data-status=inprogress],table[data-status=Finished]{display:none}"), Modernizr.localstorage && (window.sessionStorage.tennislivetabs = "3");
                }
                cssblock(ls.o.tennis.style2, "#ls-tennis-list .hdr{display:none}"), s = "";
                for (var i = ls.o.tennis.groups.length; --i >= 0; ) {
                    var v = $(ls.o.tennis.groups[i]).find(".sui-draw-fix:visible");
                    if (v.length === 0) {
                        s += "#ls-tennis-list #" + ls.o.tennis.groups[i].id + " .hdr {display:none!important}"
                    }
                }
                cssblock(ls.o.tennis.style2, s + " .sui-draw-fix-gp{visibility:visible!important}");
                $(ls.o.tennis.con).removeClass("ajax-loader");
                $('#ls-no-matches').show();
            });
            if (Modernizr.localstorage && window.sessionStorage.tennislivetabs) {
                ls.o.tennis.filters.find("input").eq(window.sessionStorage.tennislivetabs).trigger("click")
            } else {
                ls.o.tennis.filters.find("input:eq(1)").trigger("click")
            }
        };
        (function(p) {
            if (p && typeof p == "string") {
                tables(p);
                ls.o.trs.reverse()
            }
            if (p && typeof p == "object") {
                var i = p.length;
                while (--i >= 0) {
                    p[i] && tables(p[i])
                }
                ls.o.trs.reverse()
            }
        })(p.tables);
        p.scores && scores(p.scores);
        p.betting && betting(p.betting);
        p.game && game(p.game);
        p.tennisscores && tennisscores(p.tennisscores);
        p.golfscores && golfscores(p.golfscores);
    }, update: function(p) {
        f1tables = function() {
            function process(d) {
                var json = {};
                json.laps = d.race.laps;
                json.cars = d.race.car;
                json.session = d.race.session.sessionTitle;
                for (var i = json.cars.length; --i >= 0; ) {
                    e = json.cars[i];
                    var temp = {p: e['@rank'], desc: e.raceComment, lastlap: e.currentLapTime, pits: e.numberOfPits, time: e.timeOffset, s: e.position, f1desc: true}
                    json[e.driver['@id']] = temp;
                }
                return json;
            }
            p.tables = process(p.f1tables);
            tables();
        }
        function customSort(a, b) {
            return(parseInt(a.getElementsByTagName('td')[0].innerHTML) - parseInt(b.getElementsByTagName('td')[0].innerHTML));
        }
        tables = function() {
            ls.o.data.tables = p.tables;
            var rows = ls.o.trs;
            var i = rows.length;
            while (--i >= 0) {
                var lsrow = ls.o.data.tables[parseInt(rows[i].id.replace("t", "").split('-')[0], 10)];
                r = lsrow && lsrow.cr == rows[i].id.split('-')[1] ? lsrow : false;
                if (r && parseInt(rows[i].getAttribute("rel"), 10) !== parseInt([r.p]), 10) {
                    var td = $$t("td", rows[i]);
                    var h = parseInt(rows[i].getAttribute("h"), 10), posfromtop = h * (parseInt(rows[i].getAttribute('rel'), 10)), newpos = (h * parseInt(r.p, 10));
                    var y = (newpos - posfromtop);
                    newpos > posfromtop && (rows[i].className += " down");
                    newpos < posfromtop && (rows[i].className += " up");
                    y === 0 && (rows[i].className = rows[i].className.replace(/ up| down|up|down/gi, ''));
                    if (!!Modernizr.csstransitions) {
                        rows[i].style.cssText = "-moz-transform:translate(0," + (y) + "px);-webkit-transform:translate(0," + (y) + "px);-ms-transform:translate(0," + (y) + "px);-o-transform:translate(0," + (y) + "px);";
                        rows[i].setAttribute("y", y);
                    }
                    function cc(p, c, e) {
                        for (var i = p.length; --i >= 0; ) {
                            if (p[i].className.match(c) && e) {
                                p[i].innerHTML = e;
                                break;
                            }
                        }
                    }
                    cc(td, '-pos', r.s || r.p);
                    cc(td, '-pts', r.pts);
                    cc(td, '-pl', r.pl);
                    if (r.f1desc && r.desc != "" && td[3]) {
                        td[3].setAttribute('colspan', '3');
                        td[3].className = '-desc ixt';
                        td[4] && td[4].parentNode.removeChild(td[4]);
                        td[4] && td[4].parentNode.removeChild(td[4]);
                    }
                    cc(td, '-desc', r.desc);
                    cc(td, '-gd', r.gd);
                    cc(td, '-gs', r.gs);
                    cc(td, '-time', r.lastlap);
                    cc(td, '-diff', r.time);
                    ls.o.data.tables.laps && ($('.standings-laps').html(ls.o.data.tables.laps));
                    ls.o.data.tables.session && ($('.standings-session').html(ls.o.data.tables.session.replace(/qualifying/i, 'Qual')))
                }
                if (rows[i] === rows[i].parentNode.getElementsByTagName('tr')[0] && !Modernizr.csstransitions) {
                    var rowstr = $(rows[i].parentNode).find('tr'), tb = rowstr[0].parentNode, order = rowstr.sort(customSort), l = order.length, f = document.createDocumentFragment();
                    for (j = 0; j < l; j++) {
                        order[j] && f.appendChild(order[j]);
                    }
                    tb.appendChild(f)
                }
            }
        };
        betting = function() {
            ls.o.data.betting = p.betting.events || [];
            ls.o.data.betob = {};
            for (var i = -1; ++i < ls.o.data.betting.length; ) {
                var row = ls.o.data.betting[i];
                ls.o.data.betob[row.event_id] = {display: row.attributes.displayed, h: row.markets[0].selections[0].price, a: row.markets[0].selections[2].price, d: row.markets[0].selections[1].price, status: row.markets[0].status.toLowerCase(), markets: row.more_markets, hs: row.markets[0].selections[0].status.toLowerCase(), as: row.markets[0].selections[2].status.toLowerCase(), ds: row.markets[0].selections[1].status.toLowerCase(), hprice: [row.markets[0].selections[0].live_price_numerator, row.markets[0].selections[0].live_price_denominator], dprice: [row.markets[0].selections[1].live_price_numerator, row.markets[0].selections[1].live_price_denominator], aprice: [row.markets[0].selections[2].live_price_numerator, row.markets[0].selections[2].live_price_denominator]}
            }
            ;
            function processlink(p, num, den) {
                var lpnum = /lp_num=(?:0x)?\d+&/;
                var lpden = /lp_den=(?:0x)?\d+/;
                p = p.replace(lpnum, 'lp_num=' + num + "&");
                p = p.replace(lpden, 'lp_den=' + den + "&");
                return p;
            }
            function findmove(n, o) {
                n = processodd(n);
                o = processodd(o);
                var currodds = o.match('/') ? parseInt(o.split("/")[0]) / parseInt(o.split("/")[1]) : o;
                var newodds = n.match('/') ? parseInt(n.split("/")[0]) / parseInt(n.split("/")[1]) : n;
                var c = currodds > newodds ? 'down' : 'up';
                currodds === newodds && (c = '');
                return c;
            }
            ;
            function processodd(o) {
                o = o.toString();
                o = o.match('evens') ? '1/1' : o;
                o = o.match('100/30') ? '10/3' : o;
                if (o.length < 4) {
                    return o
                }
                ;
                o = (function(o) {
                    var n = parseInt(o.split('/')[0]) / parseInt(o.split('/')[1]);
                    return n === Math.round(n) ? n.toString() : o;
                })(o);
                return o;
            }
            ;
            for (var i = ls.o.fix.length; --i >= 0; ) {
                var rid = ls.o.fix[i].id.replace('b', '');
                if (ls.o.fix[i].className.match(/ixf-betting/) && !ls.o.fix[i].className.match(/disabled/) && ls.o.data.betob[rid]) {
                    var jrow = ls.o.data.betob[rid];
                    var odds = $$t('a', ls.o.fix[i]);
                    ls.o.fix[i].setAttribute('data-status', jrow.status);
                    if (odds) {
                        jrow.display === 'Y' ? $(ls.o.fix[i]).removeClass('disabled') : $(ls.o.fix[i]).addClass('disabled');
                        odds[0].setAttribute('data-status', jrow.hs);
                        var oddsh = $$c('ix-oo', odds[0]);
                        odds[0].href = processlink(odds[0].href, jrow.hprice[0], jrow.hprice[1]);
                        odds[0].setAttribute('data-odds-move', findmove(jrow.h, $(oddsh).text()));
                        $(oddsh).text(processodd(jrow.h));
                        var oddsd = $$c('ix-oo', odds[1]);
                        odds[1].href = processlink(odds[1].href, jrow.dprice[0], jrow.dprice[1]);
                        odds[1].setAttribute('data-status', jrow.ds);
                        odds[1].setAttribute('data-odds-move', findmove(jrow.d, $(oddsd).text()));
                        $(oddsd).text(processodd(jrow.d));
                        var oddsa = $$c('ix-oo', odds[2]);
                        odds[2].href = processlink(odds[2].href, jrow.aprice[0], jrow.aprice[1]);
                        odds[2].setAttribute('data-status', jrow.as);
                        odds[2].setAttribute('data-odds-move', findmove(jrow.a, $(oddsa).text()));
                        $(oddsa).text(processodd(jrow.a));
                        odds[3].innerHTML = '+' + ls.o.data.betob[rid].markets;
                    }
                }
            }
        };
        scores = function() {
            ls.o.data.scores = p.scores;
            if (ls.o.list) {
                var i = ls.o.fix.length;
                while (--i >= 0) {
                    var rid = ls.o.fix[i].id.replace("-clone", "").replace("m", "");
                    var r = ls.o.data.scores[rid];
                    if (r && ls.o.fix[i].className.match(/ixf/) && !ls.o.fix[i].className.match(/ixf-betting/)) {
                        var div = $$t("div", ls.o.fix[i]);
                        var para = $$t("p", ls.o.fix[i]);
                        var active = ls.o.fix[i].className.match("ko") && r.s === "ip";
                        r.s === "et" && (active = true);
                        var ht = para[2] ? !para[2].innerHTML.match(/ht/i) && r.ht : false;
                        var att = para[2] ? !para[2].innerHTML.match(/att/i) && r.att : false;
                        if (active || r.s || ht || att) {
                            var ems = $$t("em", div.length === 5 ? div[2] : div[3]), xh = "", xa = "";
                            if (ems[0]) {
                                xh = ems[0].innerHTML != r.h;
                                xa = ems[1].innerHTML != r.a
                            }
                            if (r.s.match(/\bp\b/i)) {
                                r.h = "P";
                                r.a = "P";
                            }
                            var fulltable = $$id('ls-tbl-full'), tablematch = $$c('ixm', fulltable);
                            for (var tlen = tablematch.length; --tlen >= 0; ) {
                                if (tablematch[tlen].id == "tm" + rid) {
                                    var ts = $$t("span", tablematch[tlen]), hts = 0, ats = 1;
                                    tablematch[tlen].getAttribute('data-away') && (hts = 1, ats = 0)
                                    ts[hts] && (ts[hts].innerHTML = r.h);
                                    ts[ats] && (ts[ats].innerHTML = r.a);
                                    var status = r.s;
                                    "goal" == status && (status = "L");
                                    "ip" == status && (status = "L");
                                    ts[2] && (ts[2].innerHTML = status.toUpperCase());
                                }
                            }
                            var html = "<em class='" + xh + "'>" + r.h + "</em> <em class='right " + xa + "'>" + r.a + "</em>", scoresdiv = div.length === 5 ? div[2] : div[3];
                            r.h && r.a && (scoresdiv.innerHTML = html);
                            var hg = "", ag = "";
                            if (para.length > 0) {
                                var strong = r.hSyn && r.hSyn.match(':') ? "<strong>" : "", strongc = r.hSyn && r.hSyn.match(':') ? "</strong>" : "";
                                hg = r.hSyn ? strong + r.hSyn.replace(/\~/gi, "<br>" + strong).replace(/\:/gi, strongc) : "";
                                ag = r.aSyn ? strong + r.aSyn.replace(/\~/gi, "<br>" + strong).replace(/\:/gi, strongc) : "";
                                para[0] && (para[0].innerHTML = hg);
                                para[1] && (para[1].innerHTML = ag);
                                para[2] && r.ht && ($$t("span", para[2])[0].innerHTML = "HT:" + r.ht);
                                para[2] && r.att && ($$t("span", para[2])[1].innerHTML = "Att: " + r.att);
                                para[3] && r.desc && (para[3].innerHTML = r.desc)
                            }
                        }
                        r.s = r.s.replace(/\s/g, '-');
                        ls.o.fix[i].className.match(/minus/) ? ls.o.fix[i].className = "ixf minus " + r.s : ls.o.fix[i].className = "ixf " + r.s;
                        r.s && !r.s.match(/ko|ip|\bp\b|\ba\b/i) && para[0] && (div[4].innerHTML = r.s)
                    }
                }
            }
        };
        game = function() {
            var r = ls.o.data.game = p.game;
            if (r.s == "ft" || r.s == "aet") {
                ls.o.game.status.innerHTML = r.s
            }
            r.att && (ls.o.game.att.innerHTML = "<strong>Att:</strong> " + r.att);
            r.desc && (ls.o.game.desc.innerHTML = r.desc);
            r.ht && (ls.o.game.ht.innerHTML = r.ht);
            ls.o.game.hs.innerHTML != r.h ? (ls.o.game.hs.innerHTML = r.h, ls.o.game.hs.className = "ls-match-score true") : ls.o.game.hs.className = "ls-match-score";
            ls.o.game.as.innerHTML != r.a ? (ls.o.game.as.innerHTML = r.a, ls.o.game.as.className = "ls-match-score true") : ls.o.game.as.className = "ls-match-score";
            if (r.goals) {
                var hg = "", ag = "";
                for (gs in r.goals[0]) {
                    hg += gs + " " + r.goals[0][gs] + "</br>"
                }
                for (gs in r.goals[1]) {
                    ag += gs + " " + r.goals[1][gs] + "</br>"
                }
            }
            hg && (ls.o.game.hg.innerHTML = hg), ag && (ls.o.game.ag.innerHTML = ag);
            ls.o.game.score.className = r.s;
            var i = ls.o.game.players.length;
            while (--i >= 0) {
                var id = ls.o.game.players[i].id.replace("p", "");
                if (r.players[id]) {
                    var pl = r.players[id], li = ls.o.game.players[i], card = "";
                    var c = $$c("cards", li)[0];
                    var s = $$c("subs", li)[0];
                    pl.y == 1 && pl.r == 0 && (card = "y1");
                    pl.y == 2 && (card = "y2");
                    pl.r == 1 && (card = "r1");
                    c.className = "cards " + card;
                    pl.out > 0 && (s.className = "subs out", s.innerHTML = pl.out);
                    pl["in"] > 0 && (s.className = "subs in", s.innerHTML = pl["in"])
                }
            }
        };
        golfscores = function(p) {
            ls.o.golf.data = p;
            var tr = "", players = ls.o.golf.data.event.tournament.players.player
            for (var i = -1; ++i < players.length; ) {
                var p = players[i], dataid = "p" + p["@attributes"].id, trx = $$id(dataid), trxclone = $$id(dataid + "-clone");
                var td = trx.getElementsByTagName("td"), tdx = trxclone ? trxclone.getElementsByTagName("td") : [];
                setHTML([td[0], tdx[0]], (p.totals["@attributes"].tied == 0 ? (p.totals["@attributes"].position) : "T" + p.totals["@attributes"].position) || convertStatus(p.totals["@attributes"].status));
                setHTML([td[2], tdx[2]], p.totals["@attributes"].totaltopar);
                setHTML([td[3], tdx[3]], rounds(p, "holes"));
                setHTML([td[8], tdx[8]], p.totals["@attributes"].strokes);
                var prevholes = rounds(p);
                setHTML([td[4], tdx[4]], prevholes[0]);
                setHTML([td[5], tdx[5]], prevholes[1]);
                setHTML([td[6], tdx[6]], prevholes[2]);
                setHTML([td[7], tdx[7]], prevholes[3]);
                changeorder(trx, i + 1);
                trxclone && trxclone.setAttribute("data-pos", trx.getAttribute("data-pos"));
            }
            ls.reorderTR(ls.o.golf.tbody1);
            var date = new Date(ls.o.golf.data.event['@attributes'].updatetimestamp);
            $(ls.o.golf.timestamp).html(date.toString().slice(4));
            !Modernizr.csstransitions && ls.reorderTR(ls.o.golf.tbody2);
            function convertStatus(p) {
                console.log(p);
                p = p.replace(/disqualified/gi, 'DSQ').replace(/withdrawn/gi, "WD");
                return p;
            }
            function rounds(p, holes) {
                var tds = [], lasthole = "", current = "";
                if (p.round) {
                    for (var i = -1; ++i < p.round.length; ) {
                        tds[i] = p.round[i] ? p.round[i]["@attributes"].strokes : "-";
                        if (p.round[i] && p.round[i]["@attributes"].holesplayed && p.round[i]["@attributes"].holesplayed != "") {
                            lasthole = p.round[i]["@attributes"].holesplayed
                        }
                        if (p.round[i] && p.round[i]["@attributes"].totaltopar && p.round[i]["@attributes"].totaltopar != "") {
                            current = p.round[i]["@attributes"].totaltopar
                        }
                    }
                } else {
                    tds = ["-", "-", "-", "-"]
                }
                if (holes === "holes") {
                    return lasthole
                }
                if (holes === "currentscore") {
                    return current
                }
                return tds
            }
            function changeorder(tr, pos, trs) {
                if (tr) {
                    tr.setAttribute("data-pos", pos);
                    var y = !isNaN(parseInt(pos)) ? (parseInt(pos) - 1) * (bskyb.smallscreen ? 32 : 40) - parseInt(tr.getAttribute("rel"), 10) : 0;
                    y !== 0 && (tr.style.cssText = "-webkit-transform:translate(0," + y + "px);-moz-transform:translate(0," + y + "px);-o-transform:translate(0," + y + "px)")
                }
            }
            ;
            function setHTML(p, h) {
                for (var i = p.length; --i >= 0; ) {
                    p[i] && (p[i].getElementsByTagName('span')[0] ? p[i].getElementsByTagName('span')[0].innerHTML = h : p[i].innerHTML = h)
                }
            }}
        ;
        tennisscores = function() {
            function getsetscore(p, s) {
                for (var i = p.length; --i >= 0; ) {
                    if (p[i]['@Type'] === s) {
                        return p[i]['@name'];
                        break;
                    }
                }
            }
            function process(c) {
                for (var b = c.tournaments, c = [], d = [], f = {}, a = b.length; 0 <= --a; ) {
                    c.push(b[a])
                }
                for (a = c.length; 0 <= --a; ) {
                    if (c[a].length) {
                        for (var b = c[a], e = b.length; 0 <= --e; ) {
                            d.push(b[e].matches)
                        }
                    } else {
                        d.push(c[a].matches)
                    }
                }
                for (a = d.length; 0 <= --a; ) {
                    if (d[a].length) {
                        b = d[a];
                        for (e = b.length; 0 <= --e; ) {
                            f[b[e]["id"]] = b[e]
                        }
                    } else {
                        f[d[a]["id"]] = d[a]
                    }
                }
                return f
            }
            ;
            ls.o.tennis.data = process(p.tennisscores);
            for (var i = ls.o.tennis.matches.length; --i >= 0; ) {
                var table = ls.o.tennis.matches[i], id = table.id.replace('t', ''), data = ls.o.tennis.data[id];
                if (!data) {
                    continue
                }
                var sets = {s1: $$c('-set1', table), s2: $$c('-set2', table), s3: $$c('-set3', table), s4: $$c('-set4', table), s5: $$c('-set5', table), pts: $$c('-pts', table)}
                table.setAttribute('data-status', data['statusTxt']);
                if (data.comment || data['statusTxt'].match(/To Finish|Delayed|Delay|Retired/ig)) {
                    var sublinks = $(table).find('.sublinks')[0];
                    var comm = data.comment !== undefined ? data.comment : data['statusTxt'].match(/Retired/ig) ? data['statusTxt'].replace(/player/ig, 'player ') : "Match delayed";
                    if (sublinks) {
                        sublinks.innerHTML = comm;
                    }
                    else {
                        var thead = table.getElementsByTagName('thead')[0];
                        $(thead).after("<tfoot><tr><td colspan='7'><ul class='sublinks'><li>" + (comm) + "</li></ul></td></tr></tfoot>");
                    }
                }
                else {
                    $(table).find('tfoot').remove();
                }
                if (data['statusTxt'].match(/retired/gi)) {
                    table.setAttribute('data-status', 'Finished');
                }
                ;
                if (data['statusTxt'].match(/Cancelled/gi)) {
                    continue;
                }
                ;
                if (data['statusTxt'].match(/NSY/gi)) {
                    sets.pts[0].innerHTML = "Scheduled";
                }
                ;
                if (data['statusTxt'].match(/inprogress|finished|retired/gi)) {
                    data["firstToServe"] && table.setAttribute('data-server', 'server' + data["firstToServe"]);
                    function getSNS(p) {
                        if (!p) {
                            return p === 0 ? '0' : '-'
                        }
                        else {
                            return p
                        }
                    }
                    sets.s1[0].innerHTML = getSNS(data.sideOne.set['1']);
                    sets.s1[1].innerHTML = getSNS(data.sideTwo.set['1']);
                    sets.s2[0].innerHTML = getSNS(data.sideOne.set['2']);
                    sets.s2[1].innerHTML = getSNS(data.sideTwo.set['2']);
                    sets.s3[0].innerHTML = getSNS(data.sideOne.set['3']);
                    sets.s3[1].innerHTML = getSNS(data.sideTwo.set['3']);
                    sets.s4[0].innerHTML = getSNS(data.sideOne.set['4']);
                    sets.s4[1].innerHTML = getSNS(data.sideTwo.set['4']);
                    sets.s5[0].innerHTML = getSNS(data.sideOne.set['5']);
                    sets.s5[1].innerHTML = getSNS(data.sideTwo.set['5']);
                    if (!data['statusTxt'].match(/finished|retired/gi)) {
                        sets.pts[0].innerHTML = data.sideOne.points !== undefined ? data.sideOne.points : "";
                        sets.pts[1].innerHTML = data.sideTwo.points !== undefined ? data.sideTwo.points : "";
                    }
                    if (data.statusTxt.match(/finished|retired/ig)) {
                        var num = data.sideOne.winner ? 0 : 1;
                        var tr = $(table).find('tbody').find('tr');
                        $(tr).find('.-pts').html("");
                        $(tr[num]).addClass('sui-draw-win').find('.-pts').html('res');
                    }
                    table.setAttribute('hasloaded', true);
                }
            }
        }
        p.tables && tables();
        p.betting && betting();
        p.f1tables && f1tables();
        p.scores && scores();
        p.tennisscores && tennisscores();
        p.golfscores && golfscores(p.golfscores);
        p.game && game();
    }};
(function($) {
    $.fn.LiveArticle = function(options, f) {
        if(this.attr('data-url')) {
            var $el = this, d = {url: $el.attr('data-url'), id: '', list: 'ul', listEl: 'li', addClass: '', prefix: 'c', loadEl: '<li class="loading"></li>', showAllBtn: '#show-all', showAllUrl: '', lastUpdateText: $el.find('#last-update'), timestamp: $el.attr('data-timestamp'), stripe: false, wrapper: ['', ''], cachebuster: $el.attr('data-epoch'), interval: 0, intervalName: '', template: $el.attr('data-template'), onUpdate: function() {
                }, corrections: true};
            d.showAllUrl = $(d.showAllBtn).attr('data-url');
            var o = $.extend(d, options), list = $el.find(o.list).eq(0), loadEl = $(o.loadEl), showAllBtn = $el.find(o.showAllBtn).bind('click', showAllHandler), tint = false, template, children = list[0].getElementsByTagName(o.listEl);
            var lastID = children[0] && parseInt(children[0].id.split(o.prefix)[1]) || 0;
            var firstID = children[children.length - 1] && parseInt(children[children.length - 1].id.split(o.prefix)[1]);
        }
        else {
            var o = {};
            o.interval = 0;
            var $el = this;
        }
        function buildLi(line) {
            var li = '<li id="' + o.prefix + line.id + '" ', p = o.wrapper[0] + line.text + o.wrapper[1], time = '', clss = 'class="' + o.addClass + ' ';
            if (tint && o.stripe) {
                clss += 'tint ';
                tint = false;
            } else if (o.stripe && !tint) {
                tint = true;
            }
            if (line.line_type) {
                clss += line.line_type.replace(/ /gi, '-').toLowerCase();
            }
            if (line.minutes !== undefined) {
                time = '<em class="time">' + line.minutes + '</em>';
            } else {
                clss += ' no-time';
            }
            if (clss != 'class="') {
                li += clss + '"';
            }
            li += '>' + time + p + '</li>';
            return li;
        }
        function buildTemplate(line) {
            var li = template.replace(/{{id}}/ig, line.id).replace(/{{text}}/ig, line.text).replace(/{{time}}/ig, line.time || "").replace(/{{line_heading}}/ig, line.minutes).replace(/{{minutes}}/ig, (line.minutes || line.time || "")).replace(/{{line_type}}/ig, (line.line_type) ? line.line_type.replace(/ /gi, '-').toLowerCase() : 'no-icon');
            return li;
        }
        function updateSummary(summary) {
            $el.find('#line-summary').html(summary);
        }
        function appendData(join, html, delay) {
            var chunk = $(html);
            if (join == 'prepend') {
                setTimeout(function() {
                    var li = chunk.prependTo(list).addClass('hdn');
                    $el[0].offsetTop;
                    li.removeClass('hdn');
                    $el.children('div').attr('data-updated', 'true');
                }, delay);
            }
            else if (join == 'append') {
                chunk.appendTo(list).show();
            }
            loadEl.remove()
        }
        function updateLi(obj) {
            if (o.corrections) {
                var li = (o.template) ? buildTemplate(obj) : buildLi(obj);
                list.children('#' + o.prefix + obj.id).replaceWith(li);
            }
        }
        function updateTime(timestamp) {
            if (timestamp) {
                o.lastUpdateText.html(timestamp);
                $el.attr('data-timestamp', timestamp);
            }
        }
        function showAllHandler(e) {
            loadEl.appendTo(list);
            $el.update({join: 'append', live: false, url: o.showAllUrl});
            showAllBtn.remove();
            return false;
        }
        this.beforeUpdate = function() {
            loadEl.prependTo(list);
        };
        this.afterUpdate = function() {
            loadEl.fadeOut(600, function() {
                $(this).remove();
            })
            o.onUpdate(el);
        };
        this.update = function(opts) {
            var defaults = {join: 'prepend', live: true, url: o.url, timestamp: o.timestamp, data: false}, options = $.extend(defaults, opts), dataHandler = function(data) {
                var lines = data.line, l = lines.length, li, da = 0, html = '', i;
                if (o.stripe) {
                    if (options.join == 'prepend' && !list.find(o.listEl + ':first-child').hasClass('tint')) {
                        tint = true;
                    } else if (options.join == 'append' && !list.find(o.listEl + ':last-child').hasClass('tint')) {
                        tint = true;
                    }
                }
                if (options.live === true) {
                    for (i = l - 1; i >= 0; i--) {
                        if (lines[i].id > lastID) {
                            li = (o.template !== undefined) ? buildTemplate(lines[i]) : buildLi(lines[i]);
                            appendData(options.join, li, da * 2500);
                            da++;
                            lastID = lines[i].id;
                        } else {
                            updateLi(lines[i]);
                        }
                    }
                    updateTime(data.formattedLU || data.lastUpdated);
                    updateSummary(data.title);
                    if (o.cachebuster !== undefined) {
                        o.cachebuster = data.epochTime;
                    }
                } else {
                    for (i = 0; i < l; i++) {
                        if (lines[i].id < firstID) {
                            li = (o.template !== undefined) ? buildTemplate(lines[i]) : buildLi(lines[i]);
                            html += li;
                        }
                    }
                    firstID = lines[l - 1].id;
                    if (html !== '') {
                        appendData(options.join, html, 0);
                    }
                }
            };
            if (!options.updateItem) {
                if (o.cachebuster !== undefined && options.live === true) {
                    options.url = o.url.split('$').join(o.cachebuster);
                }
                $.getJSON(options.url, dataHandler);
            } else {
                dataHandler(options[options.updateItem]);
            }
        };
        if (o.interval > 0) {
            var self = this;
            o.intervalName = setInterval(function() {
                self.update();
            }, o.interval);
        }
        if (o.template) {
            template = $.trim($('#' + o.template).html());
        }
        return this;
    };
})(jQuery);
(function($) {
    var clicktype = Modernizr.touch ? 'touchstart' : 'click';
    function update(oo, l) {
        oo.lis.length > 1 && !oo.con.getAttribute('height') && (oo.con.style.height = "auto", oo.con.style.overflow = "visible");
        l && ($('#photo0').remove(), oo.railmask.style.height = oo.liheight + "px", oo.lis = $(oo.rail).find('.l'), oo.imgs = $(oo.rail).find('.img'));
        l && oo.lis.length > 1 && $(oo.con).removeClass('rail-hdn');
        oo.hash && oo.con.className.match(/rail-giant/) && hash(oo);
        oo.con.setAttribute('loaded', 'done');
        if (oo.clones) {
            oo.lis = oo.clones;
            oo.imgs = $(oo.clones).find('.img')
        }
        ;
        oo.height && (oo.con.style.height = oo.lis[oo.pos].getAttribute('height') + "px");
        oo.rail.style.width = oo.con.className.match(/t2/) ? oo.width + "px" : (oo.lis.length * oo.liwidth + "px");
        oo.steps = Math.ceil(oo.lis.length / Math.floor(oo.viewable));
        oo.count.innerHTML = (oo.pos + 1) + " <span>/</span> " + oo.steps;
        oo.imgs[oo.pos] && (oo.caption.innerHTML = oo.imgs[oo.pos].getAttribute('alt') || oo.imgs[oo.pos].getAttribute('title') || "");
        loadimage(oo, oo.pos);
        if (oo.viewable < oo.lis.length && !oo.con.className.match(/t2/)) {
            $(oo.nav).show()
        }
        if (oo.lis.length > 1) {
            if (oo.lis.length * oo.liwidth < oo.width) {
                $(oo.back).add(oo.forward).hide();
            }
            else {
                $(oo.back).add(oo.forward).show()
            }
            if (!oo.cyclical) {
                $(oo.forward).add(oo.back).removeClass('disabled');
                oo.pos === oo.steps - 1 && (oo.forward.className += ' disabled', oo.pos = oo.steps - 1);
                oo.pos === 0 && (oo.back.className += ' disabled', oo.pos = 0);
            }
            if (oo.junior && oo.page === oo.steps) {
                $(oo.forward).add(oo.back).removeClass('disabled');
                $(oo.forward).addClass('disabled');
            }
            if (oo.junior && oo.page === 1) {
                $(oo.forward).add(oo.back).removeClass('disabled');
                $(oo.back).addClass('disabled');
            }
        }
    }
    function hash(oo, p) {
        var loaded = oo.con.getAttribute('loaded');
        if (loaded && !oo.overlay) {
            location.hash = "!/" + oo.hash + "/" + oo.pos
        }
    }
    function twin(oo) {
        if (oo.twin && !oo.junior) {
            $('#' + oo.twin).find('.l').eq(oo.pos).trigger('click', 'xrail');
        }
    }
    function liactions(e, t, oo, x) {
        var tempos = oo.pos;
        var index = $(oo.lis).index(t);
        if (oo.overlay && !bskyb.smallscreen) {
            overlay(oo, index);
        }
        else {
            $(oo.lis).removeClass('on');
            t.className += ' on';
            oo.pos = index;
            if ((!oo.junior && x === 'x') || (!oo.junior && x === 'xhash')) {
                move(oo.rail, (oo.pos * oo.liwidth))
                update(oo);
            }
            if ((oo.junior && x !== 'x' && oo.twin) || (oo.junior && oo.twin && x === 'xhash')) {
                oo.page = Math.floor(oo.pos / 9) + 1;
                var x = (oo.width - oo.offset) * (oo.page - 1);
                move(oo.rail, x);
                var twin = $$id(oo.twin);
                if (twin.getAttribute('data-cyclical')) {
                    $(twin).find('.forward').trigger(clicktype, [oo.pos, 'once']);
                }
                else {
                    $(twin).find('.l').eq(oo.pos).trigger('click', 'x');
                    update(oo);
                }
                oo.pos = Math.max(oo.page - 1, 0);
            }
        }
        ;
    }
    function tools(e, oo) {
        if (e.target.id === 'rail-thumbs') {
            var s = document.getElementById(oo.twin), t = e.target;
            t.className.match(/\bon\b/) ? ($(s).removeClass('on'), t.className = "enable-tooltip") : ($(s).addClass('on'), t.className = "enable-tooltip on");
        }
        if (e.target.id === 'rail-fullscreen' && Modernizr.fullscreen) {
            var t = e.target;
            if (t.className.match(/\bon\b/)) {
                e.target.className = '';
                bskyb.browser.match(/chrome/) && document.webkitCancelFullScreen();
                bskyb.browser.match(/firefox/) && document.mozCancelFullScreen();
            }
            else {
                t.className = 'on';
                var el = $$id('rail-overlay') || $$id('content');
                bskyb.browser.match(/chrome/) && el.webkitRequestFullScreen(Element.ALLOW_KEYBOARD_INPUT);
                bskyb.browser.match(/firefox/) && el.mozRequestFullScreen();
            }
            document.addEventListener("mozfullscreenchange", function() {
                !document.mozFullScreen && ($$id('rail-fullscreen').className = '');
            }, false);
            document.addEventListener("webkitfullscreenchange", function() {
                !document.webkitIsFullScreen && ($$id('rail-fullscreen').className = '');
            }, false);
        }
        if (e.target.id === 'rail-fullscreen' && e.target.className.match(/\bon\b/)) {
        }
    }
    function navsetup(oo, f) {
        var n = document.createElement('div');
        n.className = 'rail-nav', n.style.cssText = 'display:none;clear:both', o = {}, t = "";
        var h = oo.con.offsetHeight;
        h > 197 && (t = " t2");
        h > 600 && (t = " t3");
        n.innerHTML = "<div class='sui-arw back " + t + "'><img src='/img/sui.gif' class='sui' /></div><div class='count'></div><div class='sui-arw forward " + t + "'><img src='/img/sui.gif' class='sui' /></div>";
        !$$c('caption', oo.con)[0] && (n.innerHTML += '<div class="caption"></div>');
        oo.railcon.appendChild(n);
        oo.nav = n;
        oo.back = $$c('back', n)[0];
        oo.forward = $$c('forward', n)[0];
        oo.caption = $$c('caption', oo.con)[0];
        oo.count = $$c('count', n)[0];
        f && f(oo);
    }
    function navactions(e, t, oo, f) {
        if (t.className.match(/disabled/)) {
            return false
        }
        if (t.className.match(/back/)) {
            var x = (oo.width - oo.offset) * (oo.pos - 1);
            move(oo.rail, x);
            oo.page--;
            oo.pos--;
        }
        else {
            var x = (oo.width - oo.offset) * (oo.pos + 1);
            move(oo.rail, x);
            oo.page++;
            oo.pos++;
        }
        twin(oo);
        f && f(oo);
    }
    var navcycle = 0;
    function navactionscycle(e, t, oo, x, f) {
        if (navcycle === 0) {
            navcycle = 1;
            this.back = function(t) {
                !t && oo.pos--;
                oo.pos === -1 && (oo.pos = oo.clones.length - 1);
                var n = oo.clones[oo.pos].cloneNode(true);
                oo.rail.insertBefore(n, oo.rail.firstChild);
                oo.rail.style.marginLeft = (oo.liwidth) * -1 + "px";
                $(oo.rail).animate({marginLeft: 0}, {queue: false, duration: 550, complete: function() {
                        navcycle = 0;
                        oo.rail.removeChild(n.nextSibling);
                    }});
            }
            this.forward = function(t) {
                !t && oo.pos++;
                oo.pos === oo.steps && (oo.pos = 0);
                oo.rail.appendChild(oo.clones[oo.pos]);
                $(oo.rail).animate({marginLeft: (oo.liwidth) * -1}, {queue: false, duration: 550, complete: function() {
                        navcycle = 0;
                        var l = oo.rail.getElementsByTagName('li')[0];
                        oo.rail.removeChild(l);
                        oo.rail.style.marginLeft = 0;
                    }});
            }
            if (x || x === 0) {
                if (x < oo.pos) {
                    oo.pos = x;
                    this.back('t');
                }
                else {
                    oo.pos = x;
                    this.forward('t');
                }
            }
            else {
                if (t.className.match(/back/)) {
                    this.back();
                }
                else {
                    this.forward();
                }
            }
            update(oo);
            twin(oo);
            f && f(oo);
        }
    }
    function touchswipe(oo, func) {
        var twin = $(document.getElementById(oo.twin));
        var speed = 333.33;
        $(oo.rail).swipe({triggerOnTouchEnd: true, swipeStatus: swipeStatus, threshold: 60, allowPageScroll: 'vertical'});
        function prev() {
            oo.pos = Math.max(oo.pos - 1, 0);
            scroll((oo.width - oo.offset) * oo.pos, speed);
            loadimage(oo, oo.pos);
            oo.viewable < 2 && twin.find('.l').eq(oo.pos).trigger('click', 'xhash');
            update(oo);
        }
        function next() {
            oo.pos = Math.min(oo.pos + 1, oo.steps - 1);
            scroll((oo.width - oo.offset) * oo.pos, speed);
            loadimage(oo, oo.pos);
            oo.viewable < 2 && twin.find('.l').eq(oo.pos).trigger('click', 'xhash');
            update(oo);
        }
        function swipeStatus(c, b, a, d) {
            b === "move" && (a === "left" || a === "right") ? a === "left" ? scroll((oo.width - oo.offset) * oo.pos + d, 0) : a === "right" && scroll((oo.width - oo.offset) * oo.pos - d, 0) : b === "cancel" ? scroll((oo.width - oo.offset) * oo.pos, speed) : b === "end" && (a === "right" ? prev() : a === "left" && next())
        }
        function scroll(c, b) {
            if (!oo.con.className.match(/rail-tiny/)) {
                oo.rail.style.cssText += ";-webkit-transition-duration:" + (b / 1E3).toFixed(1) + "s;";
                oo.rail.style.cssText += ";-o-transition-duration:" + (b / 1E3).toFixed(1) + "s;";
                var a = (c < 0 ? "" : "-") + Math.abs(c).toString();
                oo.rail.style.cssText += ";-webkit-transform:translate3d(" + a + "px,0,0);";
                oo.rail.style.cssText += ";-o-transform:translate(" + a + "px,0);";
            }
        }}
    function close(oo) {
        var o = $$id('rail-overlay');
        o.parentNode.removeChild(o);
        location.hash = "!";
    }
    function overlay(oo, index, f) {
        var o = document.createElement('div');
        o.id = 'rail-overlay';
        document.body.appendChild(o);
        document.body.offsetTop;
        o.className = "on";
        var gallery = document.createElement('div'), rail1 = document.createElement('div'), rail2 = document.createElement('div'), lis = "", lis2 = "";
        gallery.id = 'rail-giant-wrap';
        for (var len = oo.lis.length, i = -1; ++i < len; ) {
            lis += "<li class='l'><img src='" + (oo.imgs[i].src.replace('160', '75x75')) + "' class='img' alt='" + oo.imgs[i].alt + "'></li>";
            lis2 += "<li class='l'><img src='" + (oo.imgs[i].src.replace('160', '800x600').replace('75x75', '800x600')) + "' class='img' alt='" + oo.imgs[i].alt + "'></li>";
        }
        rail1.id = 'o-small';
        rail1.className = 'rail rail-small enable-rail';
        rail1.setAttribute('data-twin', 'o-big');
        rail1.setAttribute('data-junior', 'true');
        rail1.setAttribute('data-hash', oo.con.id);
        rail1.setAttribute('data-first-run', index);
        rail1.innerHTML = '<div class="rail-con"><div class="rail-mask"><ul class="u"> ' + lis + '</ul></div></div>';
        rail2.id = 'o-big';
        rail2.className = 'rail rail-giant enable-rail"';
        rail2.setAttribute('data-twin', 'o-small');
        rail2.setAttribute('data-hash', oo.con.id);
        rail2.setAttribute('data-first-run', index);
        rail2.innerHTML = '<div class="rail-con"><div class="rail-mask"><ul class="u"> ' + lis2 + '</ul></div></div><h1>' + oo.title + '</h1><div class="caption"> </div><div id="rail-head"><div id="rail-close" class="rail-close hdn"><img src="/img/sui.gif" id="rail-close-tile" /></div><div id="rail-thumbs" class="enable-tooltip" title="toggle thumbnails"> </div><div id="rail-fullscreen"> </div></div>';
        gallery.appendChild(rail1);
        gallery.appendChild(rail2);
        o.appendChild(gallery);
        $(rail1).rail();
        $(rail2).rail();
    }
    function loadimage(oo, pos) {
        function loadme(img) {
            if (img && img.getAttribute('data-src')) {
                var newsrc = img.getAttribute('data-src');
                img.onload = function() {
                    this.removeAttribute('data-src');
                    this.src.match(/800x600/) && callomniture(5000);
                };
                img.naturalWidth > 0 && (img.src = newsrc, img.removeAttribute('data-src'));
                img.readyState === "complete" && (img.src = newsrc, img.removeAttribute('data-src'));
                img.src = newsrc;
            }
        }
        if (oo.viewable >= 2) {
            for (var i = ((Math.floor(oo.viewable) * 2) * (pos) - 1); ++i < (Math.floor(oo.viewable) * 2) * (pos + 1); ) {
                loadme(oo.imgs[i]);
            }
        }
        else {
            oo.imgs[oo.pos - 1] && loadme(oo.imgs[oo.pos - 1]);
            oo.imgs[oo.pos + 1] && loadme(oo.imgs[oo.pos + 1]);
            loadme(oo.imgs[oo.pos]);
        }
    }
    function move(oo, x) {
        if (Modernizr.csstransitions) {
            oo.style.cssText += ";-webkit-transform:translateX(" + (x * -1) + "px);-moz-transform:translateX(" + (x * -1) + "px);-o-transform:translateX(" + (x * -1) + "px)";
        }
        else {
            $(oo).animate({"marginLeft": (x * -1) + "px"}, {duration: 500, queue: false})
        }
    }
    $.fn.rail = function(options) {
        return this.each(function(options) {
            var $this = this, oo = {};
            oo.mobile = bskyb.smallscreen;
            oo.con = $this;
            oo.railcon = $$c('rail-con', this)[0];
            oo.railmask = $$c('rail-mask', this)[0];
            oo.rail = $$c('u', this)[0];
            oo.lis = $(oo.con).find('.l');
            oo.liheight = oo.lis[0].offsetHeight;
            oo.width = oo.con.offsetWidth;
            oo.con.setAttribute('original-width', oo.width);
            oo.imgs = $(oo.con).find('.img');
            oo.smallscreen && oo.con.className.match(/rail-giant/) && (function() {
                for (var i = -1; ++i < oo.imgs.length; ) {
                    var d = oo.imgs[0].getAttribute('data-src');
                    d = d.replace(/800x600/, '330');
                    oo.imgs[0].setAttribute('data-src', d);
                }
            })();
            oo.height = oo.con.getAttribute('data-height') || false;
            (function() {
                for (var i = oo.lis.length; --i >= 0; ) {
                    oo.height && oo.lis[i].setAttribute('height', oo.lis[i].offsetHeight);
                }
            })();
            oo.con.className.match(/rail-tiny/) && oo.lis.length > 1 && (oo.lis[0].style.visibility = 'visible');
            oo.caption = $$c('caption', this);
            oo.count = $$c('count', this);
            oo.id = oo.con.id;
            oo.tools = $(oo.con).find('#rail-head');
            oo.liwidth = oo.lis[0].offsetWidth;
            oo.viewable = oo.width / oo.liwidth;
            oo.steps = Math.ceil(oo.lis.length / Math.floor(oo.viewable));
            oo.page = 1;
            oo.offset = oo.width % oo.liwidth;
            oo.railmask.style.width = (oo.liwidth * oo.viewable) + "px";
            oo.hash = oo.con.getAttribute('data-hash') || false;
            oo.firstrun = oo.con.getAttribute('data-first-run') || false;
            if (oo.hash && !oo.firstrun) {
                oo.firstrun = location.hash.split('/')[2]
            }
            ;
            oo.touch = Modernizr.touch;
            oo.pos = parseInt($this.getAttribute('data-first-run'), 10) || 0;
            oo.cyclical = oo.con.getAttribute('data-cyclical') || false;
            oo.touch && (oo.cyclical = false);
            oo.overlay = oo.con.getAttribute('data-overlay') || false;
            oo.twin = oo.con.getAttribute('data-twin') || false;
            oo.junior = oo.con.getAttribute('data-junior') || false;
            oo.auto = parseInt(oo.con.getAttribute('data-auto'), 10) || false;
            oo.title = oo.con.getAttribute('data-title') || false;
            navsetup(oo, function(n) {
                oo = n
            });
            if (!oo.con.id.match(/vid-more/)) {
                if (oo.con.getAttribute('data-overlay') && bskyb.smallscreen === true) {
                }
                else {
                    $(oo.con).on('click', '.l:not(.on)', function(e, x) {
                        liactions(e, this, oo, x)
                    });
                }
            }
            $(oo.con).on('click', '#rail-head', function(e) {
                tools(e, oo)
            });
            $(oo.con).on('click', '#rail-close', close);
            oo.con.className.match(/rail-giant/) && $(document).on('keyup', function(e) {
                if ($$id('o-big')) {
                    e.keyCode === 27 && close(oo);
                    e.keyCode === 37 && $(oo.back).trigger(clicktype);
                    e.keyCode === 39 && $(oo.forward).trigger(clicktype);
                }
            });
            function orientation() {
                if (navigator.userAgent.toLowerCase().match(/iphone|ipod|ipad/)) {
                    return((window.orientation / 2) % 2) === 0 ? 'portrait' : 'landscape';
                }
                else {
                    return screen.width > screen.height ? 'landscape' : 'portrait';
                }
            }
            ;
            if (bskyb.smallscreen && oo.con.className.match(/rail-story|rail-big|rail-giant/)) {
                oo.devicestyle = document.createElement('style');
                oo.devicestyle.id = "orientation";
                oo.devicestyle.setAttribute('rel', orientation());
                document.head.appendChild(oo.devicestyle);
            }
            var orienttype = navigator.userAgent.toLowerCase().match(/iphone|ipod|ipad/) ? 'orientationchange' : 'resize';
            oo.con.className.match(/rail-story|rail-big|rail-giant/) && bskyb.smallscreen && (window.addEventListener(orienttype, function() {
                oo.devicestyle.getAttribute('rel') != orientation() && (function() {
                    oo.width = document.getElementById('container').offsetWidth - 24;
                    oo.liwidth = oo.width;
                    oo.devicestyle.innerHTML = ".rail-giant .rail-mask,.rail-giant,.rail-giant .l,.rail-big .l,.rail-big,.rail-big .rail-mask,.rail-story .l,.rail.rail-story,.rail-story .rail-mask{width:" + (oo.width) + "px!important}";
                    oo.devicestyle.setAttribute('rel', orientation());
                })();
                update(oo);
            }, false))
            update(oo);
            if (oo.firstrun && oo.con.id === 'o-big') {
                $(oo.con).find('.l').eq(oo.firstrun).trigger('click', 'x');
            }
            if (oo.cyclical && oo.viewable < 2 && !oo.touch) {
                oo.clones = $(oo.lis).clone(true);
                oo.rail.innerHTML = "";
                oo.rail.appendChild(oo.clones[0]);
                $(oo.con).on(clicktype, '.sui-arw', function(e, x) {
                    navactionscycle(e, this, oo, x)
                });
            }
            else {
                oo.cyclical = false;
                oo.con.removeAttribute('data-cyclical');
                oo.touch && !oo.con.className.match(/t2/) && !oo.con.className.match(/story/) && touchswipe(oo);
                oo.touch && oo.con.className.match(/story/) && oo.lis.length > 1 && (!navigator.userAgent.match(/iPhone|iPod/i)) && (!navigator.userAgent.match(/chrome/i)) && touchswipe(oo);
                $(oo.con).on(clicktype, '.sui-arw', function(e) {
                    navactions(e, this, oo, function() {
                        update(oo);
                    })
                });
            }
            try {
                oo.rail.addEventListener('DOMNodeInserted', function() {
                    !oo.cyclical && update(oo, 'live');
                });
            } catch (e) {
                oo.railcon.onpropertychange = function(e) {
                    !oo.cyclical && update(oo, 'live');
                };
            }
            $this.auto = function() {
                if (oo.cyclical && oo.lis.length > 1) {
                    oo.interval = setInterval(function() {
                        !window.dropdownmenuisopen && $(oo.forward).trigger(clicktype)
                    }, oo.auto)
                }
                else {
                }
            };
            oo.auto && !oo.touch && $this.auto();
            window.onhashchange = function(e) {
                var hashindex = parseInt(location.hash.split('photo/')[1], 10);
                oo.hash && $('#' + oo.twin).find('.l').eq(hashindex).trigger('click', 'xrail');
            };
            (function() {
                if (oo.auto) {
                    var t = document.getElementById(oo.twin);
                    $(oo.con).add(t).mouseenter(function() {
                        clearInterval(oo.interval);
                    });
                    $(oo.con).add(t).mouseleave(function() {
                        $this.auto();
                    });
                }
            })();
        })
    }
})(jQuery)
$($$c('enable-rail')).rail();
(function() {
    if (!Modernizr.touch) {
        var d = document.getElementsByTagName("body")[0];
        $(d).delegate(".enable-tooltip", "mouseenter mousemove mouseleave", function(a) {
            var b = document.getElementsByTagName("tooltip");
            if (a.type == "mouseenter") {
                var c = document.createElement("tooltip");
                c.id = "tooltip";
                c.className = this.getAttribute("data-class");
                var text = this.getAttribute("data-text") || this.title || this.alt || false;
                if (text) {
                    if (this.getAttribute('data-url')) {
                        $(c).load(this.getAttribute('data-url'))
                    }
                    else {
                        c.innerHTML = text
                    }
                    ;
                    d.appendChild(c);
                    d.offsetTop;
                    $(b[0]).addClass('on');
                    this.getAttribute('title') && this.setAttribute("oldtitle", this.getAttribute('title'));
                    this.getAttribute('alt') && this.setAttribute("oldalt", this.getAttribute('alt'));
                    this.removeAttribute("title");
                    this.removeAttribute("alt");
                    var x = a.pageX;
                    var y = a.pageY;
                    y = (document.documentElement.clientHeight - y) + 5;
                    b[0].style.cssText = "bottom:" + y + "px;left:" + a.pageX + "px;";
                }
            }
            if (a.type == "mousemove") {
                var ww = window.innerWidth;
                var w = b[0] && b[0].offsetWidth;
                var x = a.pageX;
                var y = a.pageY;
                y = (document.documentElement.clientHeight - y) + 5;
                if ((w + x) - ww > -40) {
                    (x = x - (w + 30));
                    $(b[0]).addClass('tooltip-east');
                }
                else {
                    $(b[0]).removeClass('tooltip-east');
                }
                b[0] && (b[0].style.cssText = "bottom:" + y + "px;left:" + x + "px;");
            }
            if (a.type == "mouseleave") {
                this.getAttribute('oldtitle') && this.setAttribute('title', this.getAttribute("oldtitle"));
                this.getAttribute('oldalt') && this.setAttribute('alt', this.getAttribute("oldalt"));
                this.removeAttribute("oldtitle");
                this.removeAttribute("oldalt");
                b[0] && b[0].parentNode.removeChild(b[0])
            }
        })
    }
})();
function parseDate(a) {
    var d = /^(\d{1,2})[\.,/](\d{1,2})([\.,/](\d{1,4})?)?$/;
    try {
        var e = a.replace(/^\s+|\s+$/g, "");
        if (e.match(d)) {
            var g = e.replace(d, "$1").stripPrecedingZeroes(), h = e.replace(d, "$2").stripPrecedingZeroes(), f;
            var i = (e.replace(d, "$4") || (new Date).getFullYear().toString()).stripPrecedingZeroes(), b = parseInt(i);
            f = b > 100 ? b.toString() : b >= 70 ? "19" + b.toString() : b >= 10 ? "20" + b.toString() : "200" + b.toString();
            var c = new Date;
            c.setFullYear(parseInt(f));
            c.setMonth(parseInt(h) - 1);
            c.setDate(parseInt(g));
            if (!(c.getFullYear().toString() == f && (c.getMonth() + 1).toString() == h && c.getDate() == g))
                throw Error();
            return c
        }
    } catch (j) {
        return null
    }
}
String.prototype.toDate = function() {
    return parseDate(this)
};
String.prototype.stripPrecedingZeroes = function() {
    for (var a = 0; this.charAt(a) == "0"; )
        a++;
    return this.substring(a)
};
function isDate(a) {
    return!!parseDate(a)
}
Date.parseDate = parseDate;
String.prototype.isDate = function() {
    return isDate(this)
};
var TINY = {};
function T$(f) {
    return document.getElementById(f)
}
function T$$(f, g) {
    return g.getElementsByTagName(f);
}
TINY.table = function() {
    function f(d, func) {
        this.n = d;
        this.pagesize = 20;
        this.paginate = 0;
        this.func = func;
    }
    function g(d) {
        d = T$(d);
        d.b = T$$("tbody", d)[0];
        d.r = d.b.rows;
        var slidy = $('#slide-info');
        slidy.remove();
        $(d).find('a').removeClass('on');
        return d
    }
    function h(d, c) {
        var i;
        var e;
        var a, b;
        e = a = d.v.toLowerCase(), d = e;
        i = b = c.v.toLowerCase(), c = i;
        a.match(/\//) && a.match(/\//g).length > 1 ? (a = parseDate(d), b = parseDate(c)) : (a.match(/\//) && a.match(/\//g).length == 1 ? (a = parseInt(a.split("/")[0]) / parseInt(a.split("/")[1]), a = parseFloat(a)) : a = parseFloat(d.replace(/(\$|\,)/g, "")), b.match(/\//) && b.match(/\//g).length == 1 ? (b = parseInt(b.split("/")[0]) / parseInt(b.split("/")[1]), b = parseFloat(b)) : b = parseFloat(c.replace(/(\$|\,)/g, "")));
        isNaN(a) && isNaN(b) && (a = d, b = c);
        return a > b ? 1 : a < b ? -1 : 0
    }
    f.prototype.init = function(d, c) {
        var a = g(d), b = 0;
        this.parent = a;
        this.e = d;
        this.l = a.r.length;
        a.a = [];
        a.h = T$$("thead", T$(d))[0].rows[0];
        for (a.w = a.h.cells.length; b < a.w; b++) {
            var e = a.h.cells[b];
            if (!e.className.match("nosort")) {
                var img = new Image;
                img.src = '/img/sui.gif';
                img.className = 'sui sui-sort';
                e.appendChild(img);
                e.onclick = new Function(this.n + ".wk(this.cellIndex)")
            }
        }
        for (b = 0; b < this.l; b++) {
            a.a[b] = {}
        }
        c != null && (new Function(this.n + ".wk(" + c + ")"))();
        if (this.paginate) {
            this.g = 1, this.pages()
        }
    };
    f.prototype.wk = function(d) {
        for (var c = g(this.e), a = c.h.cells[d], b = 0; b < this.l; b++) {
            c.a[b].o = b;
            var e = c.r[b].cells[d];
            var hdn = $$c('hdn', e)[0];
            var sortme = $$c('sortme', e)[0];
            c.a[b].v = (hdn && hdn.innerHTML) || (sortme && sortme.innerHTML) || (e.innerText ? e.innerText : e.textContent);
            c.a[b].v = c.a[b].v.match(/\bevs/i) ? '1' : c.a[b].v;
        }
        for (b = 0; b < c.w; b++) {
            if (e = c.h.cells[b], !e.className.match("nosort")) {
            }
        }
        c.p == d ? (c.a.reverse(), c.d ? $(a).removeClass('desc').addClass('asc') : $(a).removeClass('asc').addClass('desc'), c.d = c.d ? 0 : 1) : ($(c.h).find('th').removeClass('asc').removeClass('desc'), c.p = d, c.a.sort(h), c.d = 0, $(a).removeClass('desc').addClass('asc'));
        this.th = a;
        a = document.createElement("tbody");
        for (b = 0; b < this.l; b++) {
            e = c.r[c.a[b].o].cloneNode(true);
            a.appendChild(e);
            b % 2 == 0 ? $(e).removeClass('odd').addClass('even') : $(e).removeClass('even').addClass('odd');
            for (var e = T$$("td", e), f = 0; f < c.w; f++) {
            }
        }
        c.replaceChild(a, c.b);
        this.paginate && this.size(this.pagesize);
        this.func && this.func(this.parent, this.th, d);
    };
    f.prototype.page = function(d) {
        var c = g(this.e), a = 0, b = d + parseInt(this.pagesize);
        if (this.currentid && this.limitid) {
            T$(this.currentid).innerHTML = this.g
        }
        for (; a < this.l; a++) {
            c.r[a].style.display = a >= d && a < b ? "" : "none";
        }
    };
    f.prototype.move = function(d, c) {
        var a = d == 1 ? c ? this.d : this.g + 1 : c ? 1 : this.g - 1;
        if (a <= this.d && a > 0) {
            this.g = a, this.page((a - 1) * this.pagesize)
        }
    };
    f.prototype.size = function(d) {
        this.pagesize = d;
        this.g = 1;
        this.pages();
        this.page(0);
        if (this.currentid && this.limitid) {
            T$(this.limitid).innerHTML = this.d
        }
    };
    f.prototype.pages = function() {
        this.d = Math.ceil(this.l / this.pagesize)
    };
    return{sorter: f}
}();
(function() {
    var b = $$c("enable-sort"), c = b.length;
    for (sortabletables = {}; --c >= 0; ) {
        var def = b[c].getAttribute('data-defaultsort');
        def === 'false' ? def = null : def = b[c].getAttribute('data-sortindex') ? parseInt(b[c].getAttribute('data-sortindex')) : 0;
        var d = b[c].id, a = d.replace("-", "") + "sorter";
        sortabletables[a] = new TINY.table.sorter("sortabletables['" + a + "']", function(a, b, c) {
            c == 0 ? b.className.match(/desc/g) ? $(a).addClass("nohr") : $(a).removeClass("nohr") : $(a).addClass("nohr")
        });
        sortabletables[a].init(d, def)
    }
})();
(function($) {
    $.fn.morePlease = function(options) {
        var defaults = {list: "#" + $(this).attr("id"), handler: $(this).next("a"), basket: parseInt($(this).attr("data-basket")) || 1, maxAjaxRequests: parseInt($(this).attr("data-maxrequests")) || 9999, nextPage: parseInt($(this).attr("data-currentpage")) || 0, url: $(this).attr("data-url") || "/", locationHash: "page", locationHashHistory: true, alertMessage: $(this).attr("data-alertmessage") || "<strong>No more articles.</strong>You have reached the end of this article listing.", alertCallToAction: $(this).attr("data-alertaction") || "", alertClass: $(this).attr("data-alertclass") || "sui-msg"};
        var options = $.extend(defaults, options);
        var handler = options.handler;
        var list = options.list;
        var maxAjaxRequests = options.maxAjaxRequests;
        var nextPage = options.nextPage;
        var alertClass = options.alertClass;
        var alertMessage = options.alertMessage;
        var alertCallToAction = options.alertCallToAction;
        var url = options.url;
        var basket = options.basket;
        var locationHash = options.locationHash;
        var pagingCount = parseInt($(list).children().length);
        $(document).ready(function() {
            var pathname = window.location.hash;
            var fullLocationHash = '#!' + locationHash + '=';
            if (pathname.indexOf(fullLocationHash) > -1) {
                var ExitPageOffset = pathname.replace(fullLocationHash, "");
            }
        });
        $(handler).click(function(e) {
            e.preventDefault();
            if (maxAjaxRequests == nextPage) {
                $(list).after("<div class='" + alertClass + "'><p>" + alertMessage + "</p>" + alertCallToAction + "</div>");
                $(handler).remove();
                return;
            }
            nextPage = nextPage + 1
            $.ajax({cache: true, url: url, data: {page: nextPage, basketId: basket}, dataType: "html", beforeSend: function() {
                    $(handler).addClass("ajax-loader");
                    return;
                }, success: function(html) {
                    if ($(html).find(list).length > 0) {
                        var returnCount = $(html).find(list).children().length;
                        $(list).append($(html).find(list).html());
                        $(list).attr("data-currentPage", nextPage);
                        url = $(html).find(list).attr("data-url");
                        console.log($(html).find(list).attr("data-url"));
                        $(list).attr("data-url", $(html).find(list).attr("data-url"));
                        $(handler).removeClass("ajax-loader");
                        if (returnCount < pagingCount) {
                            ShowEndMessage();
                        }
                    } else {
                        ShowEndMessage();
                    }
                    function ShowEndMessage() {
                        $(handler).remove();
                        $(list).after("<div class='" + alertClass + "'><p>" + alertMessage + "</p>" + alertCallToAction + "</div>");
                    }
                    if (options.locationHashHistory) {
                        $(list).children().delegate('a', 'click', function(e) {
                            e.preventDefault();
                            window.location = '#!' + locationHash + '=' + nextPage + '';
                            window.location.href = $(this).attr("href");
                        });
                    }
                }});
        });
    };
})(jQuery);
(function() {
    function build(data, bar) {
        switch (bar.type) {
            case'duel':
                var lis = {}
                lis[data.entity[0]] = $$c('s-item', bar.ul[0]);
                lis[data.entity[1]] = $$c('s-item', bar.ul[1]);
                for (var j = data.items.length; --j >= 0; ) {
                    if (lis[data.entity[0]][j]) {
                        data1 = data.items[j][data.entity[0]];
                        data2 = data.items[j][data.entity[1]];
                        datatotal = parseInt(data1) + parseInt(data2)
                        w1 = parseInt((100 / datatotal) * parseInt(data1)), w2 = parseInt((100 / datatotal) * parseInt(data2));
                        w1 === 0 && (w1 = 1);
                        w2 === 0 && (w2 = 1);
                        $$c('val', lis[data.entity[0]][j])[0].innerHTML = data1;
                        $$c('val', lis[data.entity[1]][j])[0].innerHTML = data2;
                        $$c('b', lis[data.entity[0]][j])[0].style.width = (w1 * 1.5) + "%";
                        $$c('b', lis[data.entity[1]][j])[0].style.width = (w2 * 1.5) + "%";
                    }
                }
                break;
            case'bar':
                break;
            }
    }
    for (var stat = $$c("enable-stat"), i = stat.length; --i >= 0; ) {
        var url = stat[i].getAttribute('data-url');
        url && (function(url, stat) {
            var oo = {type: stat.getAttribute('data-type') || 'duel', url: url, con: stat, ul: $$t('ul', stat), interval: parseInt(stat.getAttribute('data-interval')) || 60000}
            var inte = setInterval(function() {
                $.ajax({type: "GET", url: oo.url, cache: false, success: function(data) {
                        if (data.refresh == 'true') {
                            build(data, oo)
                        }
                        else {
                            clearInterval(inte)
                        }
                    }})
            }, oo.interval)
        })(url, stat[i])
    }
})();
(function($) {
    $.fn.brightCove = function(options) {
        var defaults = {videoId: null, container: "video", playerId: "1225427762001", publisherId: "165012893", wmode: "transparent", autoStart: true, isVid: true, isUI: true, dynamicStreaming: true, videoSmoothing: true, jsPath: "http://admin.brightcove.com/js/", experience: null};
        var options = $.extend(defaults, options);
        var jsFile = "BrightcoveExperiences.js";
        if (options.experience == "all") {
            jsFile = "BrightcoveExperiences_all.js";
        }
        ;
        return this.each(function() {
            obj = $(this);
            $.getScript(options.jsPath + jsFile, function() {
                var object, params = {}, player = brightcove.createElement("object");
                params.width = obj.width();
                params.height = obj.height();
                params.playerID = options.playerId;
                params.publisherID = options.publisherId;
                params.wmode = options.wmode;
                params.isVid = options.isVid;
                params.isUI = options.isUI;
                params.dynamicStreaming = options.dynamicStreaming;
                params.videoSmoothing = options.videoSmoothing;
                params.autoStart = options.autoStart;
                params.videoId = options.videoId;
                player.id = options.container;
                for (var i in params) {
                    object = brightcove.createElement("param");
                    object.name = i;
                    object.value = params[i];
                    player.appendChild(object);
                }
                var playerContainer = document.getElementById(player.id);
                obj.wrap("<div class='vid-wrap'/>");
                brightcove.createExperience(player, playerContainer, true);
            });
        });
    };
})(jQuery);
function makeCanvas(d) {
    if (!!Modernizr.canvas) {
        var e = document.createElement("canvas"), f = e.getContext("2d"), a = d.width, b = d.height;
        e.width = a;
        e.height = b;
        f.drawImage(d, 0, 0);
        a = f.getImageData(0, 0, a, b);
        for (b = 0; b < a.height; b++)
            for (var g = 0; g < a.width; g++) {
                var c = 4 * b * a.width + 4 * g, h = (a.data[c] + a.data[c + 1] + a.data[c + 2]) / 3;
                a.data[c] = h;
                a.data[c + 1] = h;
                a.data[c + 2] = h
            }
        f.putImageData(a, 0, 0, 0, 0, a.width, a.height);
        d.src = e.toDataURL();
        d.style.visibility = 'visible';
    }
}
;
location.hash === "#showgrid" && $("#container").append("<img src='/img/util/grid.gif' id='grid' style='position:absolute;top:0;z-index:100000'/>");
(Modernizr.localstorage && window.localStorage.friendlyrefresh && ($('html,body').scrollTop(window.localStorage.friendlyrefresh), window.localStorage.removeItem('friendlyrefresh')));
;
function callomniture(tolerance, noreload) {
    function call() {
        if (window.SKY_TRACKING && window.sky) {
            !noreload && (SKY_TRACKING.autoReloadedContent = true);
            sky.tracking && sky.tracking.launch(SKY_TRACKING);
        }
    }
    tolerance = tolerance || 60000;
    if (!document.body.getAttribute("data-omniture-start")) {
        document.body.setAttribute("data-omniture-start", (new Date).getTime());
        call()
    }
    else {
        var old = parseInt(document.body.getAttribute("data-omniture-start"), 10);
        if ((new Date).getTime() > old + tolerance) {
            document.body.removeAttribute("data-omniture-start")
        }
    }
}
;
(function() {
    var v = $('.video-datesearch'), vix = v.find('.ix'), input = v.find('.sui-input'), valid = new RegExp(v.attr('data-regex'), 'i') || /^(0[1-9]|[12][0-9]|3[01])[- /.](0[1-9]|1[012])[- /.](20[0][9]|201[0-9])$/, url = v.attr('data-url') || '/racing/video-form/';
    v.on(Modernizr.touch ? 'touchstart' : 'click', '.sui-ico-date', function() {
        if (this.className.match(/\bon\b/)) {
            this.className = this.className.replace(/\bon\b/gi, '');
            vix[0].className = vix[0].className.replace(/\bon\b/gi, '');
        }
        else {
            input[0].focus();
            input[0].value.length === 0 && input.addClass('invalid');
            this.className += ' on';
            vix[0].className += ' on';
        }
    })
    input.on('keyup', function(e) {
        this.value.match(valid) ? $(this).removeClass('invalid') : $(this).addClass('invalid');
    })
    $('#video-datesearch').on('submit', function(e) {
        e.preventDefault();
        if (!input[0].className.match(/\binvalid\b/)) {
            location.href = 'http://' + location.hostname + url + (input[0].value.replace(/\//gi, '-')) + '/';
        }
    })
})();
$('#sfl').on("click", function(e) {
    e.preventDefault();
    var el = $('#footer-links');
    var textEl = $(this).find("a");
    if (!el.hasClass("show")) {
        el.addClass("show");
        textEl.text("less");
    } else {
        el.removeClass("show");
        textEl.text("more");
    }
    ;
});
(function() {
    var tables = document.getElementsByTagName('table');
    for (var i = tables.length; --i >= 0; ) {
        var tds = tables[i].getElementsByTagName('td');
        for (var j = tds.length; --j >= 0; ) {
            var newcolspan = tds[j].getAttribute('data-mobile-colspan');
            bskyb.smallscreen && newcolspan && tds[j].setAttribute('colspan', newcolspan);
        }
    }
})();
(function() {
    if (bskyb.smallscreen) {
        var s = $$id('search'), i = $(s).find('.ix'), input = $(s).find('input'), n = $$id('container');
        $(s).on('touchstart', 'h4', function() {
            $(n).hasClass('on') && $('#nav-mobile-open').trigger('touchend');
            if (s.className.match(/\bon\b/)) {
                i.hide();
                s.className = s.className.replace(/\bon\b/ig, '');
            }
            else {
                i.show();
                s.className += ' on';
            }
        })
    }
})();
(function() {
    if (bskyb.smallscreen && $("body#football-live").length > 0) {
        $('#tables-list').remove();
    }
})();
(function() {
    var fbLoaded;
    var twLoaded;
    var inner;
    $('#share #fb, #footer-share #f-fb').on("click", function() {
        closeWindow();
        inner = $("#fb").find('.wrapper');
        inner.addClass("show");
        if (!fbLoaded) {
            inner.append("<div id='fb-root'/><img class='sui x' alt='Close' src='/img/sui.gif' />");
            $.getScript("//connect.facebook.net/en_US/all.js#xfbml=1").done(function() {
                inner.children('.inner').removeClass("load");
                FB.init({xfbml: true});
            })
        }
        fbLoaded = true;
    });
    $('#share #tw,  #footer-share #f-tw').on("click", function() {
        closeWindow();
        inner = $("#tw").find('.wrapper');
        inner.addClass("show");
        if (!twLoaded) {
            inner.append("<img class='sui x' alt='Close' src='/img/sui.gif' />");
            $.getScript("//platform.twitter.com/widgets.js").done(function() {
                inner.children('.inner').removeClass("load");
            });
        }
        twLoaded = true;
    });
    $("#share").delegate(".x", "click", function() {
        closeWindow();
    });
    function closeWindow() {
        var el = $("#share li .wrapper");
        if (el.hasClass("show")) {
            el.removeClass("show");
            $("#container").removeClass("fade");
        }
    }
    ;
})();
function switchbutton() {
    if ($(".switch input[type='checkbox']:checked").length < 1) {
        $('.switch').trigger('click');
    }
    $('.switch').on('click', function() {
        if ($(this).hasClass('yes')) {
            $(this).find('.input_switch').animate({'margin-left': '0px'}, 150, function() {
                $(this).parent().removeClass('yes').addClass('no');
                $(this).next('span').text('No');
                $(this).find('input[type="checkbox"]').val("no");
            });
        } else if ($(this).hasClass('no')) {
            $(this).find('.input_switch').animate({'margin-left': '38px'}, 150, function() {
                $(this).parent().removeClass('no').addClass('yes');
                $(this).next('span').text('Yes');
                $(this).find('input[type="checkbox"]').val("yes");
            });
        }
    });
}
;
switchbutton();
(function() {
    var b = $$id("search"), a = $$t("input", b)[0], c = 'placeholder'in document.createElement("input"), inputs = $(a).add('.sui-search-auto');
    inputs.parents('form').on('submit', function(e) {
        e.preventDefault();
        var q = $(this).serialize().replace('search=', '');
        var url = (location.href.split(location.hostname)[0] + location.hostname + '/search/' + q + '/1/');
        q.length > 0 && (location.href = url);
    })
    var f = b && b.getAttribute("data-url");
    if (a) {
        var inputs = $(inputs).autocomplete(f, {append: true, scrollHeight: 320, submitonClick: true, selectFirst: false, max: bskyb.smallscreen || Modernizr.touch ? 4 : 10, dataType: "jsonp", formatItem: function(p) {
                return p.toString()
            }, formatResult: function(p) {
                return p.toString()
            }})
    }
    $(inputs).result(function(a, b) {
        $(this).parents('form').submit();
    });
    !c && (a.value = a.getAttribute("placeholder"), $(a).click(function() {
        this.value = ""
    }));
    $(a).focus(function() {
        $(this).removeClass("on").addClass("on");
    });
    $(a).blur(function() {
        if (this.value.length === 0) {
            a.className = a.className.replace(/\bon\b/g, ""), !c ? a.value = a.getAttribute("placeholder") : a.value = ""
        }
    });
})();
(function() {
    var b = $$c("ixform");
    if (a = b.length > 0) {
        for (; --a >= 0; ) {
            for (var c = $$t("a", b[a]), d = c.length; --d >= 0; ) {
                c[d].className += " enable-tooltip"
            }
        }
    }
})();
(function() {
    $($$id('page')).on('click', '.a-block,.disabled', function(e) {
        e.preventDefault();
        this.className.match(/a-block/) && this.href && (location.href = this.href);
    })
})();
(function() {
    var b = $$c("click");
    if (b.length > 0) {
        $(document.getElementById('page')).delegate(".click", "click", function(b) {
            !b.target.nodeName.match(/a/i) && !b.target.className.match(/ixa/) && (a = $$c("ixa", this)[0], a && (location.href = a.href))
        })
    }
})();
(function() {
    $('.enable-video span').click(function() {
        var vEl = $(this).parent("div").attr("id");
        var vId = $(this).parent("div").attr("data-vId");
        var newId = "v" + vEl.replace(/[a-zA-Z_\s]/g, '');
        var h = $(this).height();
        var ooyalaPlayerToken = "";
        $.ajax({
            url: "/ajax/video-authenticate",
            async: false,
            success:function(data) {
                if (data)
                    ooyalaPlayerToken = data;
            },
            dataType : 'text'
        });
        $(this).prev("img").remove();
        $(this).remove();
        $("#" + vEl).prepend('<div id="' + newId + '" style="height:' + h + 'px"/>');
        var videoPlayer = OO.Player.create(newId, vId, { 
            autoplay : true,
            "analytics": {
                "tags" : ["slweb"]
            },
            "google-ima-ads-manager": {
                "showInAdControlBar": true     
            },
            embedToken: ooyalaPlayerToken
        });
    });
})();
(function() {
    var ims = $$c('greyscale');
    if (ims) {
        for (var i = ims.length; --i >= 0; ) {
            makeCanvas(ims[i]);
        }
    }
})();
(function() {
    $('.sui-show').click(function(e) {
        e.preventDefault();
        var handler = $(this);
        var url = $(this).attr("data-url");
        if (url != undefined) {
            function slideDown() {
                $('html, body').animate({scrollTop: handler.offset().top - 90}, 'slow');
                handler.parent().children('.p-stat-container').slideDown();
            }
            if (!(handler).hasClass("sui-hide")) {
                if (!(handler).hasClass("loaded")) {
                    $.ajax({url: url, success: function(html) {
                            handler.parent().append(html);
                            $(handler).addClass('sui-hide loaded');
                            slideDown();
                        }})
                } else {
                    slideDown();
                    $(handler).addClass('sui-hide');
                }
            } else {
                $(this).parent().children('.p-stat-container').slideUp("slow", function() {
                    $(handler).removeClass('sui-hide');
                });
            }
        }
    });
})();
(function() {
    var ums = {};
    els = $$c("update-me");
    for (var i = els.length; 0 <= --i; ) {
        var id = els[i].id;
        ums[id] = {};
        ums[id].refresh = ums[id].refresh || (ums[id].refresh = parseInt(els[i].getAttribute("data-refresh"), 10)), ums[id].el = els[i];
        ums[id].url = ums[id].url || (ums[id].url = els[i].getAttribute("data-url"));
        if (!ums[id].url) {
            continue
        }
        if (ums[id].el.getAttribute('data-multi')) {
            ums.multi = ums[id].el.getAttribute('data-multi');
            ums.placeholder = document.createElement('div');
            ums.placeholder.style.cssText = 'position:absolute;visibility:hidden;left:-10000px';
            ums.placeholder.setAttribute('pseudo-id', id);
            document.body.appendChild(ums.placeholder);
        }
        ums[id]["int"] = setInterval(function() {
            var div = ums[id].el;
            if (ums.multi) {
                div = ums.placeholder
            }
            $(div).load(ums[id].url, function(d) {
                if (ums.multi) {
                    var m = ums.multi.split('|');
                    for (var j = m.length; --j >= 0; ) {
                        $('#' + m[j].split('=')[0]).html($(ums.placeholder).find('#' + m[j].split('=')[1]).clone());
                    }
                }
                var t = this.getElementsByTagName("*")[0], pid = this.id || this.getAttribute('pseudo-id');
                t.getAttribute("data-kill") && clearInterval(ums[pid]["int"]);
                t.getAttribute("data-url") && (ums[pid].url = t.getAttribute("data-url"));
                t.getAttribute("data-refresh") && (ums[pid].refresh = t.getAttribute("data-refresh"));
                this.style.opacity = "1";
                callomniture(ums[id].refresh);
            })
        }, ums[id].refresh)
    }
})();
(function() {
    $(document.getElementById('container')).on((Modernizr.touch ? 'touchstart' : 'click'), '.jump-to', function(e) {
        e.preventDefault();
        var j = this.getAttribute('data-to') || this.getAttribute('href'), el = $(j), o = el[0] ? el.offset().top : false;
        o < 13 && (o = 12);
        o && $('html,body').animate({scrollTop: (o - 12)}, 400);
    });
})();
(function() {
    $(document.getElementById('page-wrap')).on('click', '.enable-showhide', function(e) {
        var target = this.getAttribute('data-target'), alttext = this.getAttribute('data-alttext'), el = $(target);
        if (this.className.match('showhideon')) {
            this.innerHTML = "Show more";
            this.className = this.className.replace('showhideon', '');
            $($$c('hdn', el[0])).hide();
        }
        else {
            this.innerHTML = alttext;
            this.className += ' showhideon';
            $($$c('hdn', el[0])).show();
        }
    });
    (function() {
        var c = $$c('enable-loadall');
        for (var i = c.length; --i >= 0; ) {
            $(c[i]).on('click', function(e) {
                var c = this.parentNode.className.match(/\bon\b/) || this.className.match(/\bon\b/);
                var button = this.getAttribute('data-nav') ? [] : this;
                button.className += ' ajax-loader';
                e.preventDefault();
                if (!c) {
                    var t = $(this.getAttribute('data-target'));
                    t.addClass('ajax-loader');
                    t.load(this.getAttribute('data-url'), function() {
                        $(this).removeClass('ajax-loader');
                        $(button).remove();
                    });
                    if (this.getAttribute('data-nav')) {
                        var l = this.parentNode.parentNode.getElementsByTagName('li');
                        $(l).removeClass('on');
                        this.parentNode.className += ' on';
                    }
                }
            })
        }
    })();
    (function() {
        function findPos(a) {
            var b = 0, curtop = 0;
            if (a.offsetParent) {
                do
                    b += a.offsetLeft, curtop += a.offsetTop;
                while (a = a.offsetParent);
                return[b, curtop]
            }
        }
        var els = $$c('enable-showall');
        $(els).on('click', function() {
            var target = this.getAttribute('data-target'), parent = $(target)[0], top = findPos(parent)[1]
            parent && (this.style.display = 'none');
            $('html,body').animate({scrollTop: top - 12}, 333, function() {
            });
            parent && (parent.className += ' active');
        })
    })();
})();
(function() {
    var filters = $$c("fltr");
    $(filters).on("click", "li", function(e) {
        var p = this.parentNode;
        var target = $(p.getAttribute('data-target'));
        var fltr = this.getAttribute("data-class").replace(/fltr-/ig, ',.fltr-');
        fltr = fltr.replace(/^\,/, '');
        $(p).find("li").removeClass("on");
        $(this).addClass("on");
        if (fltr === "all") {
            target.find(".fltr-item").show()
        } else {
            target.find(".fltr-item").hide(), target.find(fltr).show()
        }
    })
})();
$("#top-links-all").on((Modernizr.touch ? 'touchstart' : 'click'), function(e) {
    var a = $$id('top-links');
    this.className.match(/\bon\b/) ? (a.className = '', this.className = '') : (a.className = 'on', this.className = 'on');
});
(function() {
    $('#cal-filter input:not("#chk-all")').on("change", function() {
        var filter = $(this).val();
        if ($(this).attr("checked")) {
            $(this).closest('li').addClass("on");
            if ($("#cal-filter input:checked").length == 1) {
                $('#cal .item').addClass('hdn');
            }
            $('#cal .item.' + filter).removeClass('hdn');
        } else {
            $(this).closest('.on').removeClass("on");
            $('#cal .item.' + filter).addClass('hdn');
            if ($("#cal-filter input:checked").length == 0) {
                $('#cal .item').removeClass('hdn');
            }
        }
    });
})();
(function() {
    $(document.getElementById('page')).on('click', '.comment-reply', function(e) {
        var el = this;
        $('#post-comment').before('<a href="#post-comment" class="comment-reply sui-btn">Reply to this post</a>');
        $("#post-comment").insertAfter(el);
        $('#comment-form textarea').val("");
        $('#comment-form .sui-input-txt').removeClass("form-ui-valid form-ui-invalid");
        $("#server-msg").remove();
        $('#comment-Id').val($(el).closest('li').attr('id') || '');
        $(el).remove();
        e.preventDefault();
    });
    $('#comment-block').on("focus", '#comment-form', function() {
        var $this = $(this);
        if ($(this).data('clicked')) {
        }
        $(this).css({"min-height": "100px"});
        $(this).find('.hdn').fadeIn();
        $(this).data('clicked', true);
        switchbutton();
        var team = ["Arsenal", "Aston Villa", "Blackburn Rovers", "Bolton Wanderers", "Chelsea", "Everton", "Fulham", "Liverpool", "Manchester City", "Manchester United", "Newcastle United", "Norwich City", "Queens Park Rangers", "Stoke City", "Sunderland", "Swansea City", "Tottenham Hotspur", "West Brom", "Wigan Athletic", "Aberdeen", "Celtic", "Dundee United", "Dunfermline Athletic", "Hearts", "Hibernian", "Inverness CT", "Kilmarnock", "Motherwell", "Rangers", "St Johnstone", "St Mirren", "Barnsley", "Birmingham City", "Blackpool", "Brighton", "Bristol City", "Burnley", "Cardiff", "Coventry", "Crystal Palace", "Derby", "Doncaster", "Hull", "Ipswich", "Leeds", "Leicester", "Middlesbrough", "Millwall", "Nottingham Forest", "Peterborough", "Portsmouth", "Reading", "Southampton", "Watford", "West Ham"];
        var sport = ["Football", "Racing", "Cricket", "Rugby Union", "Rugby League", "Formula One", "Greyhounds", "Golf", "Tennis", "Boxing", "Snooker", "NFL", "Darts"];
        $('#fav-team').autocomplete(team, {append: true, matchContains: true, resultsClass: "ac_results"});
        $('#fav-sport').autocomplete(sport, {append: true, matchContains: true, resultsClass: "ac_results"});
    });
    $('#comment-block').on('submit', function(e) {
        var $this = $(this).find('#comment-form');
        e.preventDefault();
        console.log(1);
        $.ajax({type: 'POST', url: $this.attr("action"), data: $this.serialize(), dataType: 'json', success: function(data) {
                if (data.errors) {
                    $("#server-msg").remove();
                    $this.children('ul').find("li:first").prepend('<div id="server-msg" class="error"><h4>Your Form Has Errors:</h4><ul></ul></div>');
                    $.each(data.errors, function(i, item) {
                        $('#server-msg ul').append("<li><img src='/img/misc/warning.png' width='18px'/>" + item.error.msg + "</li>");
                        $('[name="' + item.error.name + '"]').removeClass("form-ui-valid").addClass("form-ui-invalid");
                    });
                }
                if (data.success) {
                    $this.children('ul').find("li:first").prepend('<div id="server-msg" class="success"><h4>' + data.success.msg + '</h4></div>');
                    $('#comment-form textarea, #comment-form .sui-input-txt').val("");
                    $('#comment-form .sui-input-txt').removeClass("form-ui-valid form-ui-invalid");
                    $('.validity-alert-wrapper').hide();
                    $('#comment-form textarea').css({"min-height": "55px"}).data('clicked', false);
                    $this.find('.hdn').hide();
                }
            }});
    });
})();
(function() {
    var b = document.body, friendlyrefresh, r = b.getAttribute('data-refresh');
    if (r && !bskyb.smallscreen) {
        var r = parseInt(r, 10);
        $.idleTimer(r);
        $(document).bind("idle.idleTimer", function() {
            friendlyrefresh = setTimeout(function() {
                Modernizr.localstorage && (window.localStorage.friendlyrefresh = $(document).scrollTop());
                location.href = location.href;
            }, 2000)
        });
        $(document).bind("active.idleTimer", function() {
            clearTimeout(friendlyrefresh);
        });
    }
})();
$(".enable-tbl-slide").on("click", ".enable-tbl-slideinfo", function(e) {
    e.preventDefault();
    var row = this.parentNode.parentNode, table = row.parentNode.parentNode, slide = $$id("slide-info"), $this = this, url = this.getAttribute("data-url");
    function remove(p) {
        p.parentNode.parentNode.parentNode.removeChild(p.parentNode.parentNode)
    }
    if ($this.className.match(/\bon\b/)) {
        $this.className = $this.className.replace(/\bon\b/, "");
        $(slide).slideUp(166, function() {
            remove(this)
        })
    } else {
        function newtr() {
            var newtr = document.createElement("tr");
            newtr.className = "tbl-slide-info";
            newtr.innerHTML = "<td colspan='" + row.getElementsByTagName("td").length + "' class='ixt'><div id='slide-info' class='ajax-loader'></div></td>";
            row.parentNode.insertBefore(newtr, row.nextSibling);
            url && $("#slide-info").load(url);
            $this.className += " on";
            $("#slide-info").slideDown(166)
        }
        slide ? $(slide).slideUp(166, function() {
            var thisrow = this.parentNode.parentNode.previousSibling;
            thisrow.className = thisrow.className.replace(/\bnobordercollapse\b/, "");
            remove(this);
            $(table).find(".on").removeClass("on");
            newtr()
        }) : newtr()
    }
});
(function() {
    var filters = $$c('filters-radio');
    $(filters).on("click", "input", function(e) {
        $(filters).find("li").removeClass("on");
        $(this).parents("li").addClass("on");
    })
})();
(function() {
    var printcard = $$id('print-card'), tools = $$id('racecard-tools'), racecardvid = $$id('racecard-video'), racecardaudio = $$id('racecard-audio');
    if (!Modernizr.flash) {
        if (racecardvid && !racecardvid.className.match(/disabled/)) {
            racecardvid.className += ' disabled';
            racecardvid.innerHTML = racecardvid.getElementsByTagName('a')[0].innerHTML;
        }
        if (racecardaudio && !racecardaudio.className.match(/disabled/)) {
            racecardaudio.className += ' disabled';
            racecardaudio.innerHTML = racecardaudio.getElementsByTagName('a')[0].innerHTML;
        }
    }
    $(printcard).on('click', function() {
        window.print();
    });
    tools && (tools.className += ' on');
})();
(function() {
    var horse1 = $$id("horse1-option"), horse2 = $$id("horse2-option"), sel1 = $$id("horse1"), sel2 = $$id("horse2"), form = $$id("horse-h2h"), button = $$t("button")[0], con = $$id("ajax-con");
    function enable(p) {
        var h2h = p.parentNode.getAttribute("data-h2h");
        if (h2h) {
            $(horse2).find("ul").load(h2h, function(data) {
                if (data.match("error")) {
                    sel2.className += " disabled";
                    $(sel2).find("span").html("No matches")
                } else {
                    $(sel2).find("span").html("Select a horse")
                }
            });
            if (h2h) {
                var els = $$t("a", horse2);
                sel2.className = sel2.className.replace(/disabled/, "");
                $(sel2).find("span").html("Select a horse");
                var head = $("span", sel2)[0];
                var value = $(horse2).find("optgroup").html();
                button.disabled = true;
                head.innerHTML = value;
                for (var i = els.length; --i >= 0; ) {
                    var m = els[i];
                    var h = m.getAttribute("href");
                    h.match(h2h) ? m.parentNode.className = m.parentNode.className.replace(/disabled/g, "") : m.parentNode.className += " disabled"
                }
            }
        }
    }
    $(horse1).on("click", "a", function(e) {
        enable(this)
    });
    $(horse2).delegate("a", "click", function(e) {
        !this.parentNode.className.match(/disabled/) && (button.disabled = false, button.innerHTML = "<span class='more'>Compare</span>");
        $('#h2h-submit').prop('disabled', false);
        enable(this)
    });
    $(form).submit(function(e) {
        button.className += ' ajax-loader';
        e.preventDefault();
        con.className += " ajax-loader";
        var horse1Id = $('input[name="horse1"]', this).val();
        var horse2Id = $('input[name="horse2"]', this).val();
        $(con).load("/ajax/racing/headtohead/stats/" + horse1Id + "/" + horse2Id + "", function() {
            con.style.opacity = "1";
            $(button).removeClass('ajax-loader');
            con.className = con.className.replace(/\bajax-loader\b/g)
        })
    });
})()