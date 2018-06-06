dojo.registerModulePath("dojo.cookie", "dojo/cookie");
dojo.require("dojo.cookie");

function isIOSDevice() {
    return ( navigator.userAgent.match(/iPhone/i) || navigator.userAgent.match(/iPod/i) || navigator.userAgent.match(/iPad/i) );
}

function isIPadDevice() {
    return navigator.userAgent.match(/iPad/i);
}

function isAndroidDevice() {
    return navigator.userAgent.match(/Android/i);
}

function isBlackBerryDevice() {
    return ( navigator.userAgent.match(/BlackBerry/i) || navigator.userAgent.match(/PlayBook/i) );
}

function isWindowsPhone() {
    return navigator.userAgent.match(/Windows Phone/i);
}

/*
Displays a predefined "you can dload the iPad app" message.
message - the message to display
link - the link to open if user clicks OK (normally an Apple-store link).
*/
function showIpadAppMessage( message, link, dontShowAfterNumberOfVisits, showAgainAfterNumberOfDays ) {
    if (!message || !link || !isIPadDevice()){
        return;
    }

    // get the cookie value, which is how many times we already showed the message
    var cookieName = 'ipadCookie';
    numberOfTimesPopupShown = dojo.cookie(cookieName);

    var isVisit = true;
    if (typeof numberOfTimesPopupShown === "undefined" || isNaN(numberOfTimesPopupShown)) {
        // cookie not set or expired
        numberOfTimesPopupShown = 0;
    } else {
        // cookie already set
        numberOfTimesPopupShown = Number(numberOfTimesPopupShown);

        // now look if we come from our webpage itself (referrer). We don't won't to bother the user at every click
        // rather than every visit
        if(typeof document.referrer !== "undefined" && document.referrer != "" && document.referrer.indexOf(window.location.host) != -1) {
            var isVisit = false;
        }
    }

    if(isVisit && numberOfTimesPopupShown < dontShowAfterNumberOfVisits) {
        dojo.cookie('ipadCookie', numberOfTimesPopupShown + 1, {
            expires: showAgainAfterNumberOfDays,
            path: '/'
        });

        doDownload = window.confirm(message.replace("\'","'"));
        if ( doDownload ) {
            window.location.href=link;
        }
    }
}

/*
Hides the default no-flash message.
divId = the DIV id of the msg
*/
function hideIOSFlashMessage(divId){
    if (isIOSDevice() || isAndroidDevice() || isBlackBerryDevice() || isWindowsPhone()){
        document.getElementById(divId).style.display = "none";
    }
}
