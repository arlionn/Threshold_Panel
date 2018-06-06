/**
 * @package OMS Widget : Right Box A/RightBox-Widget
 * @author Tony Simbine
 * @Version 2009.03.26.02.TS10, updated 2011.02.16 16:42
 * Copyright (c) 2009-2011 and developed by Tony Simbine <ejantos (at) yahoo.com>
 * Copyright (c) 2009-2011 OMS Online Marketing Service GmbH & Co. KG
 * @link http://www.oms.eu/
 **/
if (typeof bc_id == 'undefined') var bc_id = new Array();

function bc_widget_rightboxa(arg) {
    var id = bc_id.length;
    bc_id.push(new bc_make_widget_rightboxa(arg, id));
    bc_id[id].generate();
    bc_id[id].retrieve();
}

function bc_make_widget_rightboxa(arg, id) {
    this.id = id;
    this.postfix = arg.postfix || '';
    this.token = arg.token || (this.id > 0? bc_id[0].token: 'myReadToken.');
    this.url = arg.url || '/videocenter-Url';
    this.dir = arg.dir || './bcw_rightboxa';
    this.playerId = (!arg.playerIds) ? null : arg.playerIds[0];
    this.playlistId = (!arg.playlistIds) ? null : arg.playlistIds[0];
    this.width = arg.width || 170;
    this.xmlUrl = arg.xml || false;
    this.xml = '';
    this.httpRequest = null;
    this.videoCount = null;
    this.page = 0;
    this.sliding = false;
    this.referName = 'rightboxa';
    this.ie = (document.all) ? true: false;
    this.TStab = new Array(null);
    if (arg.TsTabs)	{
        if (arg.TsTabs[0] != null && arg.TsTabs[0].length > 1)
            this.TStab[0] = arg.TsTabs[0];
    }
    this.TStarget = null;
    this.maximum = 11;
    if (arg.TSopt)	{
        if (arg.TSopt['target'] && arg.TSopt['target'].length > 2)	{
            this.TStarget = ' target="' + arg.TSopt.target + '"';
            this.TStargets = arg.TSopt.target;
        }
        if (arg.TSopt['maximum'] && arg.TSopt['maximum'] > 0)
            this.maximum = arg.TSopt['maximum'];
    }
    this.url += ((this.url.indexOf('?') > -1)? '&': '?');
    this.tsSize = (arg.tsSize && Number(1/arg.tsSize)) || 0.75;

    this.generate = function() {
        var code = document.createElement('link');
        code.setAttribute('href', this.dir + this.postfix + '/brightcove.css');
        code.setAttribute('rel', 'stylesheet');
        code.setAttribute('type', 'text/css');
        document.getElementsByTagName('head')[0].appendChild(code);

        if(this.ie) {
            document.write('<style type="text/css">');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_title' + this.postfix + ' { height: 23px; }');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_mainVideo_title' + this.postfix + ' { height: 48px; }');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_mainVideo_thumb' + this.postfix + ' div.bc_overlay' + this.postfix + ' { background: none; filter: progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + this.dir + '/inactive_overlay_large.png",sizingMethod="image"); }');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_mainVideo_thumb' + this.postfix + ' div.bc_mainHover' + this.postfix + ' { background: none; filter: progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + this.dir + '/overlay_large.png",sizingMethod="image"); }');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_arrow_left' + this.postfix + ' { top: 29px; }');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_arrow_right' + this.postfix + ' { top: 9px; }');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_outer' + this.postfix + ' div.bc_inner' + this.postfix + ' div.bc_video' + this.postfix + ' div.bc_thumb' + this.postfix + ' { width: 84px; height: 64px; }');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_outer' + this.postfix + ' div.bc_inner' + this.postfix + ' div.bc_video' + this.postfix + ' div.bc_thumb' + this.postfix + ' div.bc_smallHover' + this.postfix + ' { background: none; filter: progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + this.dir + '/overlay_small.png",sizingMethod="scale"); }');
            document.write('div.bc_widget_rba' + this.postfix + ' div.bc_outer' + this.postfix + ' div.bc_inner' + this.postfix + ' div.bc_video' + this.postfix + ' div.bc_video_title' + this.postfix + ' { height: 40px; }');
            document.write('</style>');
        }
        document.write('<div id="' + this.id + 'BCWRA_widget" class="bc_widget_rba' + this.postfix + '">');
        document.write('<div class="bc_title' + this.postfix + '" id="' + this.id + 'BCWRA_title"></div>');
        document.write('<div id="' + this.id + 'BCWRA_thumb" class="bc_mainVideo_thumb' + this.postfix + '" onmouseover="bc_id[' + this.id + '].mainHover(1)" onmouseout="bc_id[' + this.id + '].mainHover(0)"><img id="' + this.id + 'BCWRA_mainSrc" class="bc_src' + this.postfix + '" /><div class="bc_overlay' + this.postfix + '" id="' + this.id + 'BCWRA_mainOverlay"></div></div><div class="bc_mainVideo_title' + this.postfix + '" id="' + this.id + 'BCWRA_mainTitle"></div>');
        document.write('<div id="' + this.id + 'BCWRA_arrowLeft" class="bc_arrow_left' + this.postfix + '" onclick="bc_id[' + this.id + '].slide(0)"><img src="' + this.dir + this.postfix + '/arrow_left.png" onmouseover="bc_id[' + this.id + '].arrowOver(this,0)" onmouseout="bc_id[' + this.id + '].arrowOff(this,0)" /></div><div id="' + this.id + 'BCWRA_arrowRight" class="bc_arrow_right' + this.postfix + '" onclick="bc_id[' + this.id + '].slide(1)"><img src="' + this.dir + this.postfix + '/arrow_right.png" onmouseover="bc_id[' + this.id + '].arrowOver(this,1)" onmouseout="bc_id[' + this.id + '].arrowOff(this,1)" /></div><div class="bc_outer' + this.postfix + '" id="' + this.id + 'BCWRA_outer"><div class="bc_inner' + this.postfix + '" id="' + this.id + 'BCWRA_inner">');
        document.write('</div></div>');
        document.write('<div class="bc_footer' + this.postfix + '"><a href="' + this.url + 'refer=' + this.referName + '" id="' + this.id + 'BC_TS_All"' + (this.TStarget? this.TStarget: '') + '>Alle Videos &raquo;</a></div>');
        document.write('</div>');
    }

    this.retrieve = function() {
        if(this.xmlUrl)
            this.manual();
        else
            this.auto();
    }

    this.auto = function() {
        var code = document.createElement('script');
        code.setAttribute('src', 'http://api.brightcove.com/services/library?command=find_playlist_by_id&playlist_id=' + this.playlistId + '&token=' + this.token + '&playlist_fields=id,name,videos&video_fields=id,name,thumbnailURL,length,publishedDate,videoStillURL&callback=bc_id[' + this.id + '].autoReturn');
        code.setAttribute('charset', 'UTF-8');
        code.setAttribute('type', 'text/javascript');
        document.getElementsByTagName('head')[0].appendChild(code);
    }

    this.autoReturn = function(pData) {
        document.getElementById(this.id + 'BCWRA_title').innerHTML = '<a href="' + this.url +  'bcpid=' + this.playerId + '&bclid=' + this.playlistId + '&refer=' + this.referName + '"' + (this.TStarget? this.TStarget: '') + '>' + (this.TStab[0] || pData.name) + '</a>';
        var el = document.getElementById(this.id + 'BCWRA_inner');
        var videos = pData['videos'];
        var maximum = this.maximum;
        if(videos.length < maximum)
            maximum = videos.length;
        this.videoCount = maximum - 1;
        document.getElementById(this.id + 'BCWRA_mainSrc').setAttribute('src', videos[0]['videoStillURL']);
        document.getElementById(this.id + 'BCWRA_mainTitle').innerHTML = videos[0]['name'];
        var a = this.url;
        var b = this.playerId;
        var c = this.playlistId;
        var d = this.referName;
        document.getElementById(this.id + 'BCWRA_thumb').onclick = function() {
            window.location = a + 'bcpid=' + b + '&bclid=' + c + '&bctid=' + videos[0].id + '&refer=' + d;
        }
        for(var i = 1; i < maximum; i++) {
            var code = document.createElement('div');
            code.className = 'bc_video' + this.postfix;
            code.innerHTML = '<div class="bc_thumb' + this.postfix + '" onclick="window.location=\'' + this.url + 'bcpid=' + this.playerId + '&bclid=' + this.playlistId + '&bctid=' + videos[i].id + '&refer=' + this.referName + '\';" onmouseover="bc_id[' + this.id + '].smallHover(' + i + ',1)" onmouseout="bc_id[' + this.id + '].smallHover(' + i + ',0)"><img class="bc_src' + this.postfix + '" src="' + videos[i]['thumbnailURL'] + '" /><div class="bc_overlay' + this.postfix + '" id="' + this.id + 'BCWRA_smallOverlay' + i + '"></div></div><div class="bc_video_title' + this.postfix + '">' + videos[i]['name'] + '</div>';
            el.appendChild(code);
        }
        if (this.playerId && this.playlistId)
            document.getElementById(this.id + 'BC_TS_All').href = this.url + 'bcpid=' + this.playerId + '&bclid=' + this.playlistId + '&refer=' + this.referName;
        this.count();
        this.position();
        this.arrows();
    }

    this.manual = function(pCount) {
        if(window.XMLHttpRequest) {
            this.httpRequest = new XMLHttpRequest();
            if(this.httpRequest.overrideMimeType)
                this.httpRequest.overrideMimeType('text/xml');
        } else if(window.ActiveXObject) {
            try {
                this.httpRequest = new ActiveXObject('Msxml2.XMLHTTP');
            } catch (e) {
                try {
                    this.httpRequest = new ActiveXObject('Microsoft.XMLHTTP');
                }
                catch (e) {}
            }
        }
        if(!this.httpRequest)
            return false;
        var thispointer = this;
        this.httpRequest.onreadystatechange = function() {
            thispointer.handler();
        };
        this.httpRequest.open('GET', this.xmlUrl, true);
        this.httpRequest.send('');
    }

    this.handler = function() {
        var tabs, tab, videos;
        var videoIds = '';
        if (this.httpRequest.readyState == 4) {
            if (this.httpRequest.status == 200) {
                this.xml = this.httpRequest.responseXML;
                tabs = this.xml.getElementsByTagName('lineup');
                tab = tabs[0];
                document.getElementById(this.id + 'BCWRA_title').innerHTML = '<a href="' + this.url + '"' + (this.TStarget? this.TStarget: '') + '>' + (this.TStab[0] || tab.getAttribute('name')) + '</a>';
                videos = tab.getElementsByTagName('video');
                for (var i = 0; i < videos.length; i++) {
                    var video = videos[i];
                    var title = video.getElementsByTagName('title_id')[0].childNodes[0].nodeValue;
                    videoIds += title;
                    if(i + 1 < videos.length)
                        videoIds += ',';
                }
                var code = document.createElement('script');
                code.setAttribute('src', 'http://api.brightcove.com/services/library?command=find_videos_by_ids&token=' + this.token + '&video_ids=' + videoIds + '&video_fields=id,name,thumbnailURL,length,publishedDate,videoStillURL&callback=bc_id[' + this.id + '].manualReturn');
                code.setAttribute('charset', 'UTF-8');
                code.setAttribute('type', 'text/javascript');
                document.getElementsByTagName('head')[0].appendChild(code);
            }
        }
    }

    this.manualReturn = function(pData) {
        var el = document.getElementById(this.id + 'BCWRA_inner');
        var videos = pData['items'];
        var maximum = this.maximum;
        var tabs = this.xml.getElementsByTagName('lineup');
        var tab = tabs[0];
        var xml_videos = tab.getElementsByTagName('video');
        for (var i = 0; i < videos.length; i++) {
            if (videos[i] == null)	{
                videos.splice(i, 1);
                var removednode = xml_videos[i];
                removednode.parentNode.removeChild(removednode);
                i--;
            }
            if (i >= maximum)
                break;
        }
        if(videos.length < maximum)
            maximum = videos.length;
        this.videoCount = maximum - 1;
        document.getElementById(this.id + 'BCWRA_mainSrc').setAttribute('src', videos[0]['videoStillURL']);
        document.getElementById(this.id + 'BCWRA_mainTitle').innerHTML = videos[0]['name'];
        var a = this.url;
        var b = xml_videos[0].getElementsByTagName('player_id')[0].childNodes[0].nodeValue;
        var c = xml_videos[0].getElementsByTagName('lineup_id')[0].childNodes[0].nodeValue;
        var d = this.referName;
        document.getElementById(this.id + 'BCWRA_thumb').onclick = function() {
            window.location = a + 'bcpid=' + b + '&bclid=' + c + '&bctid=' + videos[0].id + '&refer=' + d;
        }
        var _playerId, __playerId, _playlistId, __playlistId;
        _playerId = b;
        __playerId = _playerId;
        _playlistId = c;
        __playlistId = _playlistId;
        for(var i = 1; i < maximum; i++) {
            var playerId = xml_videos[i].getElementsByTagName('player_id')[0].childNodes[0].nodeValue;
            var playlistId = xml_videos[i].getElementsByTagName('lineup_id')[0].childNodes[0].nodeValue;
            if (_playerId != playerId && __playerId == _playerId)
                __playerId = playerId;
            if (_playlistId != playlistId && __playlistId == _playlistId)
                __playlistId = playlistId;
            var code = document.createElement('div');
            code.className = 'bc_video' + this.postfix;
            code.innerHTML = '<div class="bc_thumb' + this.postfix + '" onclick="window.location=\'' + this.url + 'bcpid=' + playerId + '&bclid=' + playlistId + '&bctid=' + videos[i].id + '&refer=' + this.referName + '\';" onmouseover="bc_id[' + this.id + '].smallHover(' + i + ',1)" onmouseout="bc_id[' + this.id + '].smallHover(' + i + ',0)"><img class="bc_src' + this.postfix + '" src="' + videos[i]['thumbnailURL'] + '" /><div class="bc_overlay' + this.postfix + '" id="' + this.id + 'BCWRA_smallOverlay' + i + '"></div></div><div class="bc_video_title' + this.postfix + '">' + videos[i]['name'] + '</div>';
            el.appendChild(code);
        }
        if (_playerId && _playlistId && __playerId == _playerId && __playlistId == _playlistId)
            document.getElementById(this.id + 'BC_TS_All').href = this.url + 'bcpid=' + _playerId + '&bclid=' + __playlistId + '&refer=' + this.referName;
        this.count();
        this.position();
        this.arrows();
    }

    this.count = function() {
        var outer = document.getElementById(this.id + 'BCWRA_outer');
        if(this.width < 208) {
            outer.style.left = (((this.width - 20) / 2) - 43) + 'px';
            outer.style.width = 84 + 'px';
        } else if(this.width < 300) {
            this.videoCount = this.videoCount - 1;
            outer.style.left = (((this.width - 20) / 2) - 90) + 'px';
            outer.style.width = 178 + 'px';
        } else {
            this.videoCount = this.videoCount - 2;
            outer.style.left = (((this.width - 20) / 2) - 136) + 'px';
            outer.style.width = 272 + 'px';
        }
    }

    this.position = function() {
        var widget = document.getElementById(this.id + 'BCWRA_widget');
        var title = document.getElementById(this.id + 'BCWRA_title');
        var thumb = document.getElementById(this.id + 'BCWRA_thumb');
        var image = document.getElementById(this.id + 'BCWRA_mainSrc');
        var overlay = document.getElementById(this.id + 'BCWRA_mainOverlay');
        var mainTitle = document.getElementById(this.id + 'BCWRA_mainTitle');
        var arrowRight = document.getElementById(this.id + 'BCWRA_arrowRight');
        var height = this.width * this.tsSize;
        widget.style.width = (this.width - 22) + 'px';
        height = (this.width - 20) * this.tsSize;
        image.style.height = (height - 2) + 'px';
        image.style.top = '1px';
        image.style.left = '1px';
        overlay.style.width = (this.width - 22) + 'px';
        overlay.style.height = height + 'px';
        arrowRight.style.left = ((this.width - 20) - 19) + 'px';
        if(this.ie) {
            title.style.width = (this.width - 21) + 'px';
            overlay.style.left = (((this.width - 24) / 2) - 92) + 'px';
            overlay.style.top = ((height / 2) - 70) + 'px';
            thumb.style.width = (this.width - 21) + 'px';
            thumb.style.height = (height + 2) + 'px';
            image.style.width = (this.width - 25) + 'px';
            mainTitle.style.width = (this.width - 21) + 'px';
        } else {
            title.style.width = (this.width - 28) + 'px';
            thumb.style.width = (this.width - 22) + 'px';
            thumb.style.height = height + 'px';
            image.style.width = (this.width - 24) + 'px';
            mainTitle.style.width = (this.width - 20) + 'px';
        }
    }

    this.mainHover = function(pState) {
        var overlay = document.getElementById(this.id + 'BCWRA_mainOverlay');
        if(pState == 1) {
            overlay.className = 'bc_overlay' + this.postfix + ' bc_mainHover' + this.postfix;
            document.getElementById(this.id + 'BCWRA_thumb').className = 'bc_mainVideo_thumb bc_mainVideo_thumb_over';
        } else {
            overlay.className = 'bc_overlay' + this.postfix;
            document.getElementById(this.id + 'BCWRA_thumb').className = 'bc_mainVideo_thumb';
        }
    }

    this.smallHover = function(pId, pState) {
        var overlay = document.getElementById(this.id + 'BCWRA_smallOverlay' + pId);
        var pEl = overlay.parentNode;
        if (pState == 1) {
            overlay.className = 'bc_overlay' + this.postfix + ' bc_smallHover' + this.postfix;
            pEl.className = 'bc_thumb bc_thumb_over';
        } else {
            overlay.className = 'bc_overlay' + this.postfix;
            pEl.className = 'bc_thumb';
        }
    }

    this.arrows = function() {
        var arrowLeft = document.getElementById(this.id + 'BCWRA_arrowLeft');
        var arrowRight = document.getElementById(this.id + 'BCWRA_arrowRight');
        if(this.page < 1)
            arrowLeft.className = 'bc_arrow_left' + this.postfix + ' bc_arrow_disabled' + this.postfix;
        else
            arrowLeft.className = 'bc_arrow_left' + this.postfix;
        if (this.page + 1 == this.videoCount)
            arrowRight.className = 'bc_arrow_right' + this.postfix + ' bc_arrow_disabled' + this.postfix;
        else
            arrowRight.className = 'bc_arrow_right' + this.postfix;
    }

    this.arrowOver = function(pEl, pSide) {
        if (pSide == 0)
            pEl.src = this.dir + this.postfix + '/arrow_left_hover.png';
        else
            pEl.src = this.dir + this.postfix + '/arrow_right_hover.png';
    }

    this.arrowOff = function(pEl, pSide) {
        if (pSide == 0)
            pEl.src = this.dir + this.postfix + '/arrow_left.png';
        else
            pEl.src = this.dir + this.postfix + '/arrow_right.png';
    }

    this.slide = function(pDir) {
        if (!this.sliding) {
            el = document.getElementById(this.id + 'BCWRA_inner');
            if (pDir == 0) {
                if (this.page > 0) {
                    this.page--;
                    this.sliding = true;
                    this.prepSlide(el, (this.page * -94));
                }
            } else {
                if (this.page + 1 < this.videoCount) {
                    this.page++;
                    this.sliding = true;
                    this.prepSlide(el, (this.page * -94));
                }
            }
            this.arrows();
        }
    }

    this.getNum = function(pNum) {
        if (pNum) {
            if(pNum.indexOf('px') > -1)
                return parseInt(pNum.substring(0, pNum.indexOf('px')));
            else
                return parseInt(pNum);
        } else {
            return 0;
        }
    }

    this.prepSlide = function(pEl, pEnd) {
        var moveStart = this.getNum(pEl.style['left']);
        var amountToMove = pEnd - moveStart;
        var timeStart = new Date().getTime();
        var timeEnd = timeStart + 500;
        this.doSlide(this.id, pEl, 'left', amountToMove, moveStart, 500, timeEnd);
    }

    this.doSlide = function(id, pEl, pType, pDistance, pStart, pTimeTotal, pTimeEnd) {
        var currentTime = new Date().getTime();
        var timeRemaining = Math.max(0, pTimeEnd - currentTime);
        var currentMove = parseInt(pDistance - (Math.pow(timeRemaining, 3) / Math.pow(pTimeTotal, 3)) * pDistance);
        pEl.style[pType] = (pStart + currentMove) + 'px';
        if(timeRemaining > 0)
            setTimeout(function () { bc_id[id].doSlide(id, pEl, pType, pDistance, pStart, pTimeTotal, pTimeEnd); }, 10);
        else
            this.sliding = false;
    }
}
