/**
 * @package OMS Widget : Single Widget
 * @author Tony Simbine
 * @Version 2009.05.06.01.TS20, updated 2012.11.26 16:42
 * Copyright (c) 2009-2012 and developed by Tony Simbine <ejantos (at) yahoo.com>
 * Copyright (c) 2009-2012 OMS Online Marketing Service GmbH & Co. KG
 * @link http://www.oms.eu/
**/

/** UPDATE NOTICE
 *
 * File contains adaptations, these must be restored/ported on updating the code! The places in the code are marked
 * using the following tags.
 * DMN_MOD_SHARELINK: The link specified in the "TSopt.link" parameter is used as a sharing link without any
 *      modifications/additions
 * DMN_MOD_CODEDIR: Allow specification of a code directory independently from the images directory, from which the html
 *      files for the popup-iframes will be taken.
 *
 */
if (typeof now == 'undefined') var now = new Date();
var bc_firstWidget = true;
if (typeof bc_id == 'undefined') var bc_id = new Array();

function bc_widget_single(arg) {
	var id = bc_id.length;
	bc_id.push(new bc_make_widget_single(arg, id));
	bc_id[id].generate();
	bc_id[id].retrieve();
	bc_id[id].listen(id);
}

function bc_make_widget_single(arg, id) {
	this.id = id;
	this.postfix = arg.postfix || '';
	this.token = arg.token || (this.id > 0? bc_id[0].token: 'MyReadToken.');
	this.dir = arg.dir || './bcw_single';
    // DMN_MOD_CODEDIR
    this.codeDir = arg.codeDir || '/mauritius/brightcove/js_common/sw';
    // DMN_MOD_CODEDIR end
	this.videoId = (!arg.videoIds) ? null : arg.videoIds[0];
	this.playerId = (!arg.playerIds) ? null : arg.playerIds[0];
	this.width = arg.width || 300;
	this.height = null;
	this.viewed = false;
	this.ff = false;
	this.saf = false;
	this.ie = false;
	this.ie7 = false
	this.ie6 = false;
	this.first = (bc_firstWidget) ? true : false;
	this.referName = 'singlewidget';
	this.TSshow = (!arg.TSshow) ? true : arg.TSshow;
	this.TSurl = this.TSsytyle = this.TSurlText = null;
	if (arg.TSopt)	{
		if (arg.TSopt['target'] && arg.TSopt['target'].length > 2)	{
			this.TStarget = ' target="' + arg.TSopt.target + '"';
			this.TStargets = arg.TSopt.target;
		}
		if (arg.TSopt['style'] && arg.TSopt['style'].length > 1)
			this.TSsytyle = arg.TSopt['style'];
		if (arg.TSopt['url'] && arg.TSopt['url'].length > 1)
			this.TSurl = arg.TSopt['url'];
		if (arg.TSopt['urltext'] && arg.TSopt['urltext'].length > 1)
			this.TSurlText = arg.TSopt['urltext'];
	}
	this.tsSize = (arg.tsSize && Number(1/arg.tsSize))? 1/arg.tsSize: 0.75;
	this.tsPSize = (arg.tsPSize && Number(1/arg.tsPSize))? 1/arg.tsPSize: 0.75;
	this.videoRefId = (!arg.videoRefIds)? null: arg.videoRefIds[0];
	this.tsChrom = (arg.tsChrom == null)? false: arg.tsChrom;
	this.tsWidth = (arg.pwidth == null)? 480: arg.pwidth;
	this.tsHeight = Math.round(this.tsWidth * this.tsPSize);
	if (! this.tsChrom)	{
		this.tsWidth += 6;
		this.tsHeight += 56;
	}
	
	this.generate = function() {
		var browser = this.detect()[1];
		var version = this.detect()[2];
		if (browser == 'firefox')
			this.ff = true;
		else if(browser == 'safari')
			this.saf = true;
		else if(browser == 'explorer') {
			this.ie = true;
			if (version < 7)
				this.ie6 = true;
			else if(version < 8)
				this.ie7 = true;
		}
		var code = document.createElement('link');
		code.setAttribute('href', this.dir + this.postfix + '/brightcove.css');
		code.setAttribute('rel', 'stylesheet');
		// code.setAttribute('async', true);
		code.setAttribute('type', 'text/css');
		document.getElementsByTagName('head')[0].appendChild(code);
		bc_firstWidget = false;
		var xcode = '';
		if (this.ie) {
			xcode += '<style type="text/css">';
			xcode += 'div.bc_widget_single' + this.postfix + ' div.bc_overlay_off' + this.postfix + ' { background: none; filter: progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + this.dir + this.postfix + '/bc_playOverlay_Off.png", sizingMethod="image"); }';
			xcode += 'div.bc_widget_single' + this.postfix + ' div.bc_overlay_on' + this.postfix + ' { background: none; filter: progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + this.dir + this.postfix + '/bc_playOverlay_On.png", sizingMethod="image"); }';
			xcode += 'div.bc_widget_single' + this.postfix + ' div.bc_metaOverlay' + this.postfix + ' { background: none; filter: progid:DXImageTransform.Microsoft.AlphaImageLoader(src="' + this.dir + this.postfix + '/bc_metaOverlay.png", sizingMethod="scale"); }';
			if(this.ie6)
				xcode += 'div.bc_widget_single_overlay' + this.postfix + ' { background: #000; filter:alpha(opacity=90); position: absolute; }';
			if (this.ie)
				xcode += 'div.bc_widget_single_overlay' + this.postfix + ' { position: absolute; top:0px; }';
			xcode += '</style>';
			document.write(xcode);
		}
		xcode = '<div id="' + this.id + 'BCWS_widget" class="bc_widget_single' + this.postfix + '" onmouseout="bc_id[' + this.id + '].overlaySwitch(0);" onmouseover="bc_id[' + this.id + '].overlaySwitch(1);"></div>';
		document.write(xcode);
	}
	
	this.retrieve = function() {
		var code = document.createElement('script');
		var now = new Date();
		var milli = now.getTime();
		if (this.videoRefId != null && this.videoRefId.length > 1)
			code.setAttribute('src', 'http://api.brightcove.com/services/library?Rtime=' + milli + '&command=find_video_by_reference_id&reference_id=' + this.videoRefId + '&token=' + this.token + '&video_fields=name,id,shortDescription,videoStillURL,length,publishedDate,videoStillURL&callback=bc_id[' + this.id + '].autoReturn');
		else
			code.setAttribute('src', 'http://api.brightcove.com/services/library?Rtime=' + milli + '&command=find_video_by_id&video_id=' + this.videoId + '&token=' + this.token + '&video_fields=name,id,shortDescription,videoStillURL,length,publishedDate,videoStillURL&callback=bc_id[' + this.id + '].autoReturn');
		code.setAttribute('charset', 'UTF-8');
		code.setAttribute('type', 'text/javascript');
		document.getElementsByTagName('head')[0].appendChild(code);
	}
	
	this.autoReturn = function(pData) {
		if (! pData)
			return(true);
		if (! this.videoId)
			this.videoId = pData.id;
		var pSec = Math.round(pData.length / 1000);			
		var pMin = Math.floor(pSec / 60);
		pSec %= 60;
		if (pSec < 10)	pSec = '0' + pSec;
		var now = new Date();
		var milli = now.getTime();
		var code = '<img src="' + pData.videoStillURL + '" class="bc_still' + this.postfix + '" onclick="bc_id[' + this.id + '].openBox();" id="' + this.id + 'BCWS_image" />';
		code += '<div id="' + this.id + 'BCWS_overlayOn" onclick="bc_id[' + this.id + '].openBox();" class="bc_overlay_on' + this.postfix + '" style="display: none;"></div>';
		code += '<div id="' + this.id + 'BCWS_overlayOff" class="bc_overlay_off' + this.postfix + '"></div>';
		code += '<div id="' + this.id + 'BCWS_metaOverlay" onclick="bc_id[' + this.id + '].openBox();" class="bc_metaOverlay' + this.postfix + '">' + pData.name + ' <span>(' + pMin + ':' + pSec + ')</span></div>';		
		code += '<div id="' + this.id + 'BCWS_overlay" class="bc_widget_single_overlay' + this.postfix + '" onclick="bc_id[' + this.id + '].closeBox();" style="display:none;">';
		code += '<div id="' + this.id + 'BCWS_box" class="bc_widget_single_lightbox' + this.postfix + '" onclick="event.cancelBubble = true; return false;">';
		code += '<div class="bc_close' + this.postfix + '"><a href="#" onclick="bc_id[' + this.id + '].closeBox();return false;">Schlie&szlig;en</a></div>';
        // DMN_MOD_CODEDIR
        code += '<iframe name="' + this.id + 'BCWS_experience" frameBorder="0" id="' + this.id + 'BCWS_experience" src="' + this.codeDir + '/bcw_single_blank.html" width="' + (this.tsWidth + 5) + '" height="' + (this.tsHeight + 5) + '"></iframe>';
        // DMN_MOD_CODEDIR END
		code += '<div class="bc_metaLightbox' + this.postfix + '" style="width: ' + (this.tsWidth - 8) + 'px;">';
		code += '<div class="bc_title' + this.postfix + '" id="ts_bc_title_' + this.id + '">' + pData.name + '</div>';
		code += '<div class="bc_desc' + this.postfix + '" id="ts_bc_desc_' + this.id + '">' + pData.shortDescription + '</div>';
		if (this.TSurl)
			code += '<div class="bc_ts_link' + this.postfix + '" style="' + (this.TSsytyle || 'text-align: right;line-height:14px;padding:0 3px;') + '"><a href="#" onclick="bc_id[' + this.id + '].closeBox();return(window.location.href=\'' + this.TSurl + '\');">' + (this.TSurlText || 'Alle Videos') + '</a></div>';
		code += '</div></div></div>';
		document.getElementById(this.id + 'BCWS_widget').innerHTML = code;
		this.setSize();
		(document.getElementById('bc_single_dummy') || document.body).appendChild(document.getElementById(this.id + 'BCWS_overlay'));
		this.tsAutoCheck();
	}
	
	this.setSize = function() {
		this.height = this.width * this.tsSize;
		document.getElementById(this.id + 'BCWS_widget').style.height = this.height + 'px';
		document.getElementById(this.id + 'BCWS_overlayOn').style.top = Math.round((this.height / 2) - 45) + 'px';
		document.getElementById(this.id + 'BCWS_overlayOn').style.left = Math.round((this.width / 2) - 60) + 'px';
		document.getElementById(this.id + 'BCWS_overlayOff').style.top = Math.round((this.height / 2) - 45) + 'px';
		document.getElementById(this.id + 'BCWS_overlayOff').style.left = Math.round((this.width / 2) - 60) + 'px';
		if (this.ie) {
			document.getElementById(this.id + 'BCWS_widget').style.width = this.width + 'px';
			document.getElementById(this.id + 'BCWS_image').style.width = (this.width - 5) + 'px';
			document.getElementById(this.id + 'BCWS_image').style.height = (this.height - 6) + 'px';
			document.getElementById(this.id + 'BCWS_metaOverlay').style.width = (this.width - 5) + 'px';
			document.getElementById(this.id + 'BCWS_metaOverlay').style.bottom = '-1px';
		} else {
			document.getElementById(this.id + 'BCWS_widget').style.width = (this.width - 6) + 'px';
			document.getElementById(this.id + 'BCWS_image').style.width = (this.width - 6) + 'px';
			document.getElementById(this.id + 'BCWS_image').style.height = this.height + 'px';
			document.getElementById(this.id + 'BCWS_metaOverlay').style.width = (this.width - 16) + 'px';
		}
	}
	
	this.listen = function(pId) {
		var prevScroll = window.onscroll;
		var prevResize = window.onresize;
		window.onscroll = function() {
			if (typeof prevScroll == 'function')
				prevScroll();
			bc_id[pId].position();
		};
		
		window.onresize = function() {
			if(typeof prevResize == 'function')
				prevResize();
			bc_id[pId].position();
		};
	}
	
	this.overlaySwitch = function(pState) {
		if (pState == 1) {
			document.getElementById(this.id + 'BCWS_overlayOn').style.display = 'block';
			document.getElementById(this.id + 'BCWS_overlayOff').style.display = 'none';
			document.getElementById(this.id + 'BCWS_widget').className = 'bc_widget_single' + this.postfix + ' bc_widget_single_over' + this.postfix;
		} else {
			document.getElementById(this.id + 'BCWS_overlayOn').style.display = 'none';
			document.getElementById(this.id + 'BCWS_overlayOff').style.display = 'block';
			document.getElementById(this.id + 'BCWS_widget').className = 'bc_widget_single' + this.postfix;
		}
	}
	
	this.openBox = function() {
		var overlay = document.getElementById(this.id + 'BCWS_overlay');
		var box = document.getElementById(this.id + 'BCWS_box');
		var embed = document.getElementById(this.id + 'BCWS_experience');
		try {
			this.flashDisplay(true);
		} catch(e) {}
		overlay.style.display = 'block';
		box.style.display = 'block';
		embed.style.width = (this.tsWidth + 5) + 'px';
		embed.style.height = (this.tsHeight + 5) + 'px';
		this.position();
		var milli = now.getTime();
		var tsLink = '';
		if (arg.TSopt && arg.TSopt['link'] && arg.TSopt['link'].length > 0 && arg.TSopt['link'].indexOf('http') == 0)	{
			var _tsLink = arg.TSopt['link'];
            //DMN_MOD_SHARELINK:
			tsLink = '&tsLink=' + encodeURIComponent(_tsLink);
            //DMN_MOD_SHARELINK END:
		}
        // DMN_MOD_CODEDIR
		window.frames[this.id + 'BCWS_experience'].location = this.codeDir + '/bcw_single_player.html?id=' + this.id + '&videoId=' + this.videoId + '&playerId=' + this.playerId + '&refer=' + this.referName + '&tsPSize=' + this.tsPSize + '&tsWidth='+ this.tsWidth + '&tsHeight=' + this.tsHeight + '&doply=' + (this.doply? 'true': 'false') + tsLink;
        // DMN_MOD_CODEDIR END
	}
	
	this.closeBox = function() {
		var overlay = document.getElementById(this.id + 'BCWS_overlay');
		var box = document.getElementById(this.id + 'BCWS_box');
		var iFrame = window.frames[this.id + 'BCWS_experience'];
		try {
			iFrame.bc_player.pause();
			iFrame.bc_player.stop();
			iFrame.tsRemovePlayer(this.id);
		} catch(e) {}
		overlay.style.display = 'none';
		box.style.display = 'none';
		try {
			this.flashDisplay(false);
		} catch(e) {}
	}
	
	this.position = function() {
		var widget = document.getElementById(this.id + 'BCWS_widget');
		var overlay = document.getElementById(this.id + 'BCWS_overlay');
		var box = document.getElementById(this.id + 'BCWS_box');
		var embed = document.getElementById(this.id + 'BCWS_experience');
		var size = this.dimensions();
		try {
			overlay.style.top = (size[2] - 3) + 'px';
			overlay.style.left = (size[3] - 3) + 'px';
			if (embed.offsetHeight > 100) {
				box.style.top = Math.round((size[1] / 2) - (box.offsetHeight / 2)) + 'px';
				box.style.left = Math.round((size[0] / 2) - (box.offsetWidth / 2)) + 'px';
			} else {
				box.style.top = Math.round(((size[1] / 2) - (box.offsetHeight / 2)) - 206) + 'px';
				box.style.left = Math.round((size[0] / 2) - (box.offsetWidth / 2)) + 'px';
			}
			if (this.ie) {
				overlay.style.width = (size[0] + 13) + 'px';
				overlay.style.height = (size[1] + 18) + 'px';
			} else {
				overlay.style.width = size[0] + 'px';
				overlay.style.height = size[1] + 'px';
			}
		} catch(e) {}
	}
	
	this.flashDisplay = function(pOff) {
		var setState = 'visible';
		var tsActiveObject = this.id + 'BCWS_experience';
		var tsObject = document.getElementsByTagName('object');
		if (pOff)
			setState = 'hidden';
		if (tsObject != null) {
			var apNum = tsObject.length;
			for (var i = 0; i < apNum; i++)
				tsObject[i].style.visibility = setState;
			if (tsActiveObject.length > 0 && pOff)
				document.getElementById(tsActiveObject).style.visibility = 'visible';
		}
	}
	
	this.dimensions = function() {
		var BC_WINDOW_W = 0, BC_WINDOW_H = 0, BC_SCROLL_H = 0, BC_SCROLL_W = 0;
		if (typeof(window.innerWidth) == 'number') {
			BC_WINDOW_W = window.innerWidth;
			BC_WINDOW_H = window.innerHeight;
			BC_SCROLL_H = window.pageYOffset;
			BC_SCROLL_W = window.pageXOffset;
		} else if (document.documentElement && (document.documentElement.clientWidth || document.documentElement.clientHeight)) {
			BC_WINDOW_W = document.documentElement.clientWidth;
			BC_WINDOW_H = document.documentElement.clientHeight;
			BC_SCROLL_H = document.documentElement.scrollTop;
			BC_SCROLL_W = document.documentElement.scrollLeft;
		} else if (document.body && (document.body.clientWidth || document.body.clientHeight)) {
			BC_WINDOW_W = document.body.clientWidth;
			BC_WINDOW_H = document.body.clientHeight;
			BC_SCROLL_H = document.body.scrollTop;
			BC_SCROLL_W = document.body.scrollLeft;
		}
		return [ BC_WINDOW_W, BC_WINDOW_H, BC_SCROLL_H, BC_SCROLL_W ];
	}

	this.detect = function() {
		var t;
		var b = new Array(
			{string: navigator.userAgent, subString: 'Chrome', identity: 'Chrome'},
			{string: navigator.userAgent, subString: 'OmniWeb', versionSearch: 'OmniWeb/', identity: 'OmniWeb'},
			{string: navigator.vendor, subString: 'Apple', identity: 'Safari', versionSearch: 'version'},
			{prop: window.opera, identity: 'Opera'},
			{string: navigator.vendor, subString: 'iCab', identity: 'iCab'},
			{string: navigator.vendor, subString: 'KDE', identity: 'Konqueror'},
			{string: navigator.userAgent, subString: 'Firefox', identity: 'Firefox'},
			{string: navigator.vendor, subString: 'Camino', identity: 'Camino'},
			{string: navigator.userAgent, subString: 'Netscape', identity: 'Netscape'},
			{string: navigator.userAgent, subString: 'MSIE', identity: 'Explorer', versionSearch: 'MSIE'},
			{string: navigator.userAgent, subString: 'Gecko', identity: 'Mozilla', versionSearch: 'rv'},
			{string: navigator.userAgent, subString: 'Mozilla', identity: 'Netscape', versionSearch: 'Mozilla'}
		);
		var o = new Array(
			{string: navigator.platform, subString: 'Win', identity: 'windows'},
			{string: navigator.platform, subString: 'Mac', identity: 'mac'},
			{string: navigator.platform, subString: 'Linux', identity: 'linux'}
		);
		var g = function(d) {
			for (var i = 0; i < d.length; i++) {
				var s = d[i].string;
				var p = d[i].prop;
				t = d[i].versionSearch || d[i].identity;
				if (s) {
					if (s.indexOf(d[i].subString) != -1)
						return d[i].identity;
				} else if(p)
					return d[i].identity;
			}
		};

		var h = function(s, t) {
			var i = s.indexOf(t);
			if (i == -1)
				return;
			return parseFloat(s.substring(i + t.length + 1));
		};

		browser = g(b) || 'unknown';
		version = h(navigator.userAgent, t) || h(navigator.appVersion, t) || 'unknown';
		browser = browser.toLowerCase();
		os = g(o) || 'unknown';
		return new Array(os, browser, version);
	}
	this.tsAutoCheck = function() {
		var query = window.location.search.substring(1);
		this.doply = false;
		if (query.indexOf('doply=true') > -1)	{
			this.doply = true;
			this.openBox();
		}
	}
}
