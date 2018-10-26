// Highlights the entry in the navigation.
$.each($('#navcon .s .m li a'), function(val) {
        var endpoint = $(this).children("span").length > 0 ? $(this).children("span").text() : $(this).text();
        if (endpoint.toLowerCase()==$('#navcon .s').attr('data-endPoint').toLowerCase()) {
            $(this).parent().addClass('on');
        }
});

//launch betslip
$('#betting-odds-fullmarket').on('click','a.oc',function(e){
	e.preventDefault();
	if(this.getAttribute('data-status')!=='s'){

		var ev = this.id.split('oc-')[1],
		num = this.getAttribute('data-oc-price-num'),
		den = this.getAttribute('data-oc-price-den');

		window.open('https://www.skybet.com/secure/rbs?action=go_gmulti&ev_oc_id='+ev+'&lp_num='+num+'&lp_den='+den,'betslip','width=988,height=492,left='+((screen.width-988)/2)+',top='+((screen.height-492)/2));

	}
})

//add show more for hidden markets
$('#betting-odds-fullmarket').on('click','.sui-btn',function(e){

	e.preventDefault();
	var target = $$c('market-more',this.parentNode.parentNode.parentNode.parentNode)[0];

	if(this.className.match(/\bon\b/)){
		this.className = this.className.replace(/\bon\b/gi,'');
		target.style.display = "none";
	}
	else {
		this.className += ' on';
		target.style.display = "table-row-group";
	}

});

var spLiveUpdates = {};
(function(){

    spLiveUpdates.cb = {
        interval : {},uTime : {},uCount : {},kill : {},kTime : 0,lag : 60,rfs : {},
        update : function(item,url,refresh,upObj) {
                    if (typeof spLiveUpdates.cb.uCount[item]=='undefined') spLiveUpdates.cb.uCount[item] = 0;
                    if (typeof spLiveUpdates.cb.kill[item]=='undefined') spLiveUpdates.cb.kill[item] = 0;
                    if (spLiveUpdates.cb.uCount[item] == 0) spLiveUpdates.cb.uTime[item] = url.match(/(\d{5,})/);
                    //spLiveUpdates.cb.uTime[item] = url.match(/(\d+)/);
                    spLiveUpdates.cb.uCount[item]++;
                    if (spLiveUpdates.cb.kTime && spLiveUpdates.cb.uTime[item] >= (spLiveUpdates.cb.kTime+spLiveUpdates.cb.lag)) {
                         // Kills off all intervals once they have passed the lag time.
                         if (!spLiveUpdates.cb.kill[item]) spLiveUpdates.cb.kill[item] = 2;
                    }
                    $.ajax({"type":"GET", "dataType":"json", "url":url, "cache":true,  "success":function(data) {
                        spLiveUpdates.cb.interval[item]&&clearTimeout(spLiveUpdates.cb.interval[item]);
                        var r = data.refresh || refresh, c = {};
                        if (typeof spLiveUpdates.cb.rfs[item]!='undefined') r = spLiveUpdates.cb.rfs[item];
                        c[item] = data;
                        c['updateItem'] = item;
                        upObj.update(c);
                        var newUrl = spLiveUpdates.cb.getUrl(r, item, url);
                        if (spLiveUpdates.cb.kill[item] < 2) {
                            r!='kill' &&  (spLiveUpdates.cb.interval[item] = setTimeout(function(){spLiveUpdates.cb.update(item,newUrl,r,upObj)},refresh));
                        }
						callomniture(60000);
                    },
                    // Logs details to the console and sets up the interval to try again.
                    "error":function(request, status, error) {
                        spLiveUpdates.cb.interval[item]&&clearTimeout(spLiveUpdates.cb.interval[item]);
                        var r = refresh;
                        var newUrl = spLiveUpdates.cb.getUrl(r, item, url);
                        refresh!='kill' &&  (spLiveUpdates.cb.interval[item] = setTimeout(function(){spLiveUpdates.cb.update(item,newUrl,r,upObj)},refresh));
                    }
                })
        },
        // Gets the next URL with cachebuster rounded.
        getUrl : function(rRate, item, url) {
            return url;
            spLiveUpdates.cb.interval[item]&&clearInterval(spLiveUpdates.cb.interval[item]);
            var r = rRate;
            spLiveUpdates.cb.uTime[item] = parseInt(spLiveUpdates.cb.uTime[item])+(r/1000);
            var cacheBuster = Math.floor(parseInt(spLiveUpdates.cb.uTime[item]) / (r/1000)) * (r/1000);
            return url.replace(/(\d{5,})/, cacheBuster);
        }
    }
	
    olympicResults = {
        update : function(data) {
            // if successful, replace the schedule-results div html
            if (data.olympic_results.html.length > 0) {
                // replace the html
                $('div#schedule-results-promo').hide().html(data.olympic_results.html).fadeIn(1000);
            }
        }
    }

    liveRace = {
    	update : function(data) {
    		if(data.f1live.html.length) {
    			$('#standings-list').parent().html(data.f1live.html)
    		}
    	}	
    }
    
    spMatchUpdate = {
        hs : {},as : {},eTime : 0,meTime : 0,lups : 0,ko : 0,hTot : 0,aTot : 0, frmtns : {},
        mstats : [1,2,3,5,6,7,9,10,15], // Inplay status to check.
        update : function(data) {
                var pc, p1, p2, gd;
                var type = ['', '', ' pen', ' og'];
                var data = data['game'];
                if (data.s.att != undefined) {
                    document.getElementById('ls-match-att').innerHTML = data.s.att !== "" ? '<strong>Att:</strong> '+data.s.att : "";
                }
                if (data.s.comment != undefined) {
                    document.getElementById('matchComment').innerHTML = data.s.comment !== "" ? ' '+data.s.comment : "";
                }
                $('#ls-match-score').removeClass('goal');
                $('#ls-home-score').removeClass('true');
                $('#ls-away-score').removeClass('true');

                if ((data.s.status == 30 || data.s.status == 1) && !this.lups) {
                     $('#ls-match-score').removeClass('ip ft ko');
                    if (data.s.status == 1) {
                         $('#ls-match-score').addClass('ip');
                    } else {
                        $('#ls-match-score').addClass('ko');
                    }
                    spMatchUpdate.getLups(data.s['@id']);
                } else if ((data.s.status == 4 || data.s.status == 12) && !spLiveUpdates.cb.kill['game']) {
                   $('#ls-match-score').removeClass('ip ko').addClass('ft');
                   document.getElementById("ls-match-status").innerHTML = "FT";
                   spLiveUpdates.cb.kill['game'] = 1;
                   spLiveUpdates.cb.kTime = spLiveUpdates.cb.uTime['game'];
                   spLiveUpdates.cb.rfs['game'] = 20000;
                } else if (data.s.status != undefined && data && data.s.lU > this.eTime && !spLiveUpdates.cb.kill['game']) {
                    if ($.inArray(parseInt(data.s.status), [1,2,3,5,6,7,9,10,15]) > -1) $('#ls-match-score').removeClass('goal ft ko').addClass('ip');
                    if (spMatchUpdate.hTot < data.s.h.score) {
                        $('#ls-match-score').addClass('goal');
                        $('#ls-home-score').addClass('true');
                    }
                    if (document.getElementById('aggHS') && data.s.h.agg != undefined) {
                        document.getElementById('aggHS').innerHTML = data.s.h.agg;
                    }
                    document.getElementById('ls-home-score').innerHTML = data.s.h.score, spMatchUpdate.hTot = data.s.h.score;
                    if (spMatchUpdate.aTot < data.s.a.score) {
                        $('#ls-match-score').addClass('goal');
                        $('#ls-away-score').addClass('true');
                    }
                    if (document.getElementById('aggAS') && data.s.a.agg != undefined) {
                        document.getElementById('aggAS').innerHTML = data.s.a.agg;
                    }
                    document.getElementById("ls-away-score").innerHTML = data.s.a.score, spMatchUpdate.aTot = data.s.a.score;
                    document.getElementById('ls-match-ht').innerHTML = data.s.h.hT != "" ? '(ht:'+data.s.h.hT+'-'+data.s.a.hT+')' : "";
                    spMatchUpdate.hs = data.s.h.syn.replace(/~/gi, "<br>");
                    document.getElementById("ls-home-goals").innerHTML = spMatchUpdate.hs;
                    spMatchUpdate.as = data.s.a.syn .replace(/~/gi, "<br>");
                    document.getElementById("ls-away-goals").innerHTML = spMatchUpdate.as;
                    var squads = $('.ls-match-squad ul');
                    if (data.e.events != undefined && spMatchUpdate.lups) {
                        $.each(data.e.events.event, function() {
                            switch (this.type) {
                                case '5':
                                    pc = squads.find('li[id="p'+this.player['@id']+'"] a .cards');
                                    pc.hasClass('y1') ? pc.removeClass('y2').addClass('y2') : pc.removeClass('r1').addClass('r1');
                                    pc.attr('title', spMatchUpdate.fet(this.time)+"'");
                                    break;
                                case '4':
                                    pc = squads.find('li[id="p'+this.player['@id']+'"] a .cards');
                                    pc.removeClass('y1').addClass('y1');
                                    pc.attr('title',spMatchUpdate.fet(this.time)+"'");
                                    break;
                                case '12':
                                    p1 = squads.find('li[id="p'+this.player['@id']+'"] a .subs')
                                    p1.removeClass('out').addClass('out').text(spMatchUpdate.fet(this.time));
                                    p2 = squads.find('li[id="p'+this.sub['@id']+'"] a .subs');
                                    p2.removeClass('in').addClass('in').text(spMatchUpdate.fet(this.time));
                                    break;
                                default:
                                    break;
                            }
                        });
                    }
                    this.eTime = data.s.lU;
                    //spMatchUpdate.setGs();
                }
        },
        getLups : function(mId) {

            $.ajax({type:"GET", dataType:"json", url:'/ajx/football/lineups/'+mId, cache:true,  success:function(data) {
                if (data.squads != undefined) {
                    var playerCount = 0;
                    data.squads.homeTeam.squad.player && (playerCount += data.squads.homeTeam.squad.player.length);
                    data.squads.awayTeam.squad.player && (playerCount += data.squads.awayTeam.squad.player.length);
                    if ((playerCount / 2) >= 11) { // only load squads if we have 2 full teams
                        document.getElementById("ls-homeTeam").innerHTML = '';
                        document.getElementById("ls-awayTeam").innerHTML = '';
                        $.each(data.squads, function() {
                            var tn = this.homeFlag == 1 ? 'homeTeam' : 'awayTeam';
                            if ($.inArray(parseInt(this.formation), spMatchUpdate.frmtns) < 0) this.formation = '442';
                            var frmn = String(this.formation.split(''));
                            frmn = frmn.replace(/,/gi, "-");
                            $('#ls-'+tn).append('<h5>Starting Formation: <span class="ls-sq-exp">'+frmn+'<img class="sui sui-formation" src="/img/sui.gif" alt="cross"></span></h5>')
                                        .append('<div class="ls-sq-form"> </div><ul id="ls-'+tn+'-sq"></ul>')
                                        .append('<h5>Subs</h5><ul id="ls-'+tn+'-subs"></ul>');
                            // Sorts the players in rank order.
                            this.squad.player.sort(function(a, b) {
                                return a.rank - b.rank;
                            })
                            $.each(this.squad.player, function() {
                                var lip = '<li class="ls-sq" id="p'+this['@id']+'"><a href="/football/team/player/'+this['@id']+'" class="sq-ix"><strong class="sqn">'+this.shirt+'.</strong> '+this.forename+' '+this.surname+' <em class="goals"></em><em class="subs"></em><span class="cards"></span></a></li>';
                                if (this.sub == 'N') {
                                    $('#ls-'+tn+'-sq').append(lip);
                                } else {
                                    $('#ls-'+tn+'-subs').append(lip);
                                }
                            });
                        });
                    }
                    spMatchUpdate.lups = 1;
                }
            }})
        },
        fet : function(time) {
            var ta = time.split('(');
            var rTime = ta[0].split(':');
            if (ta[1]) {
                rTime = rTime[0]+'('+ta[1];
            } else {
                rTime = rTime[0];
            }
            return rTime;
        },
        setGs : function() {
                for($key in spMatchUpdate.hs) {
                    document.getElementById("ls-home-goals").innerHTML = $key+' '+spMatchUpdate.hs[$key];
                }
                for($key in spMatchUpdate.as) {
                    document.getElementById("ls-home-goals").innerHTML = $key+' '+spMatchUpdate.hs[$key];
                }
        }
    }

    fbMatchBet = {
        update : function(d) {
            var d = d['matchBet'];
            if (d.event.markets != undefined) {
                var markets = (d.event.markets);
                for(var i = markets.length;--i>=0;) {
                    var id = markets[i].market_id;
                    $('#mkt-'+id).parent().attr('data-display',markets[i].attributes.displayed);
                    //loop through selections for that market
                    var selections = markets[i].selections;
                    for(var j = selections.length;--j>=0;) {
                        var sort = selections[j].sort;
//                        var el;
                        var findme = sort + '-odds';
                        //change odds
// This assumes the file is ordered H/D/A. It might not be.
                        var odd = $$id('oc-'+selections[j].selection_id);
                        if(odd){
                            var disableme = selections[j].status === 'S' || selections[j].attributes.displayed === 'N';
                            if (disableme) { odd.setAttribute('data-status','s') };
                            if(selections[j].live_price_numerator) {
                                odd.setAttribute('data-oc-price-num',selections[j].live_price_numerator);
                            };
                            if(selections[j].live_price_denominator) {
                                odd.setAttribute('data-oc-price-den',selections[j].live_price_denominator);
                            };
//                            var el = odd.getElementsByTagName('b')[0];
                            var cur;
//                            el[sort] = odd.getElementById(findme);
                            cur = odd.getElementById(findme).innerHTML;
//                            var cur = el.innerHTML;
                            cur.match('evens') && (cur = '1/1');
                            var currodds = parseInt(cur.split("/")[0])/parseInt(cur.split("/")[1]);
                            var newodds = parseInt(selections[j].price.split("/")[0])/parseInt(selections[j].price.split("/")[1]);
                            var c = currodds > newodds ? 'up' : 'down';
                            currodds === newodds && (c = '');
                            odd.setAttribute('data-odds-move',c);
                            odd.getElementById(findme).innerHTML = selections[j].price;
//                            el.innerHTML = selections[j].price;
                        }
                    }
                }
                
                window.marketgroups =  window.marketgroups || $('.mktgrp');
                for(var i = -1;++i<window.marketgroups.length;){
                    var parent = window.marketgroups[i];
                    var mkts = $(parent).find('.mkt');
                    var display = false;
                    for(var j = -1;++j<mkts.length;){
                       if (mkts[j].getAttribute('data-display')==='Y') { display=true };
                    }
                    if(!display) { 
                        $(parent).hide();
                    } else {
                        $(parent).show();
                    };
                }
            }
        }
    }

    fbMatchBetWdw = {
        update : function(d) {
            var d = d['matchBetWdw'];
            if (d.event.market.selections != undefined) {
                $.each(d.event.market.selections, function() {
                    $('#mb-'+this.selection_id).attr('href',this.url);
                    $('#mb-'+this.selection_id+' strong').text(this.live_price_numerator+'/'+this.live_price_denominator);
                });
            }

        }
    }

    rugby_match = {
        hs : {},as : {},eTime : 0,meTime : 0,lups : 0,ko : 0,hTot : 0,aTot : 0, type : {},
        update : function(data) {
                var pc, p1, p2, gd;
                var type = ['', '', ' pen', ' og'];
                var data = data['game'];
                if (data.s.att != undefined) document.getElementById('ls-match-att').innerHTML = data.s.att !== "" ? '<strong>Att:</strong> '+data.s.att : "";
                if (data.s.comment != undefined) document.getElementById('matchComment').innerHTML = data.s.comment.length ? data.s.comment : "";
                $('#ls-match-score').removeClass('goal');
                $('#ls-home-score').removeClass('true');
                $('#ls-away-score').removeClass('true');
                if ((data.s.status == 1 || data.s.status == 2) && !this.lups) {
                     $('#ls-match-score').removeClass('ip ft ko');
                    if (data.s.status == 2) {
                         $('#ls-match-score').addClass('ip');
                    } else {
                        $('#ls-match-score').addClass('ko');
                    }
                    rugby_match.getLups(data.s['@id']);
                } else if ((data.s.status == 5 || data.s.status == 7) && !spLiveUpdates.cb.kill['game']) {
                   $('#ls-match-score').removeClass('ip ko').addClass('ft');
                   document.getElementById("ls-match-status").innerHTML = "FT";
                   spLiveUpdates.cb.kill['game'] = 1;
                   spLiveUpdates.cb.kTime = spLiveUpdates.cb.uTime['game'];
                   spLiveUpdates.cb.rfs['game'] = 20000;
                } else if (data.s.status != undefined && data && data.s.lU > this.eTime && !spLiveUpdates.cb.kill['game']) {
                    if (data.s.status != parseInt(0)) $('#ls-match-score').removeClass('goal ft ko').addClass('ip');
                    if (rugby_match.hTot < data.s.h.score) {
                        $('#ls-match-score').addClass('goal');
                        $('#ls-home-score').addClass('true');
                    }
                    document.getElementById('ls-home-score').innerHTML = data.s.h.score, rugby_match.hTot = data.s.h.score;
                    if (rugby_match.aTot < data.s.a.score) {
                        $('#ls-match-score').addClass('goal');
                        $('#ls-away-score').addClass('true');
                    }
                    document.getElementById("ls-away-score").innerHTML = data.s.a.score, rugby_match.aTot = data.s.a.score;
                    document.getElementById('ls-match-ht').innerHTML = data.s.h.hT != "" ? '(ht:'+data.s.h.hT+'-'+data.s.a.hT+')' : "";

                    // Formats the synopsis strings
                    var regexArray = {'~':"<strong>", ":":"</strong>", "\\),":") ", "\\(":" (", ",":", "};
                    for (var val in regexArray) data.s.h.syn = data.s.h.syn.replace(new RegExp(val, "gi"), regexArray[val]);
                    document.getElementById("ls-home-goals").innerHTML = "<strong>"+data.s.h.syn;
                    for (var val in regexArray) data.s.a.syn = data.s.a.syn.replace(new RegExp(val, "gi"), regexArray[val]);
                    document.getElementById("ls-away-goals").innerHTML = "<strong>"+data.s.a.syn;

                    var squads = $('.ls-match-squad ul');
                    if (data.e.event != undefined && rugby_match.lups) {
                        $.each(data.e.event, function() {
                            switch (this.type) {
                                case '5':
                                    pc = squads.find('li[id="p'+this.player['@id']+'"] a .cards');
                                    pc.removeClass('y1').addClass('y1');
                                    pc.attr('title', rugby_match.fet(this.time)+"'");
                                    break;
                                case '6':
                                    pc = squads.find('li[id="p'+this.player['@id']+'"] a .cards');
                                    pc.hasClass('y1') ? pc.removeClass('y2').addClass('y2') : pc.removeClass('r1').addClass('r1');
                                    pc.attr('title', rugby_match.fet(this.time)+"'");
                                    break;
                                case '7':
                                    p1 = squads.find('li[id="p'+this.player['@id']+'"] a')
                                    if (!p1.find('em').hasClass('out')) p1.append('<em class="subs out">'+rugby_match.fet(this.time)+'</em>');
                                    p2 = squads.find('li[id="p'+this.sub['@id']+'"] a');
                                    if (!p2.find('em').hasClass('in')) p2.append('<em class="subs in">'+rugby_match.fet(this.time)+'</em>');
                                    break;
                                default:
                                    break;
                            }
                        });
                    }
                    this.eTime = data.s.lU;
                    //match.setGs();
                }
        },
        getLups : function(mId) {

            $.ajax({type:"GET", dataType:"json", url:'/ajx/rugby-'+rugby_match.type+'/lineups/'+mId, cache:true,  success:function(data) {
                if (data.squads != undefined) {
                    document.getElementById("ls-homeSquad").innerHTML = '';
                    document.getElementById("ls-awaySquad").innerHTML = '';
                    $.each(data.squads, function() {
                        var tn = this.homeFlag == 1 ? 'homeSquad' : 'awaySquad';
                        $('#ls-'+tn).append('</div><h5>Team:</h5><ul id="ls-'+tn+'-sq"></ul>')
                                    .append('<h5>Replacements</h5><ul id="ls-'+tn+'-subs"></ul>');
                        // Sorts players by shirt number ascending.
                        this.squad.player.sort(function(a,b){return a.shirt-b.shirt});
                        $.each(this.squad.player, function() {
                            var lip = '<li class="ls-sq" id="p'+this['@id']+'"><a href="/rugby-'+rugby_match.type+'/team/player/'+this['@id']+'" class="sq-ix"><!--<strong class="sqn">'+this.shirt+'.</strong>--> '+this.forename+' '+this.surname+' <em class="goals"></em><span class="cards"></span></a></li>';
                            if (this.sub == 'No') {
                                $('#ls-'+tn+'-sq').append(lip);
                            } else {
                                $('#ls-'+tn+'-subs').append(lip);
                            }
                        });
                    });
                    rugby_match.lups = 1;
                }
            }})
        },
        fet : function(time) {
            var ta = time.split('(');
            var rTime = ta[0].split(':');
            if (ta[1]) {
                rTime = rTime[0]+'('+ta[1];
            } else {
                rTime = rTime[0];
            }
            return rTime;
        },
        setGs : function() {
                for($key in rugby_match.hs) {
                    document.getElementById("ls-home-goals").innerHTML = $key+' '+rugby_match.hs[$key];
                }
                for($key in match.as) {
                    document.getElementById("ls-home-goals").innerHTML = $key+' '+rugby_match.hs[$key];
                }
        }
    }
	
})();


/**
 * ADTECH.loadAd - Asynchronously loads an ad into a specified html target.
 *               - Requires adtech_asyncif.html to be placed on the web server.
 *
 * @param a ID value of the container where ad will be displayed.
 * @param c Actual width of the creative.
 * @param d Actual height of the creative.
 * @param b ADTECH addyn tag URL.
 * 
 */
if(window.adgroupid==undefined){window.adgroupid=Math.round(Math.random()*1000);}
var ADTECH = ADTECH || {};
ADTECH.loadAd = function (a, c, d, b) {
    a = document.getElementById(a);
    if (a != null && typeof ADTECH_ASYNCIF_PATH != "undefined") {        
        a.style.width = c + "px";
        a.style.height = d + "px";
        a.innerHTML = "";
        b = ADTECH_ASYNCIF_PATH + "#" + escape(b);
        a.innerHTML = '<iframe width="' + c + '" height="' + d + '" src ="' + b + '" frameborder="0" marginwidth="0" marginheight="0" allowtransparency="true" scrolling="no"></iframe>';    
    }
};

/* REVENUE SCIENCE AD TAG CODE */
var rsi_segs = [];
var segs_beg=document.cookie.indexOf('rsi_segs=');
if(segs_beg>=0){
    segs_beg=document.cookie.indexOf('=',segs_beg)+1;
    if(segs_beg>0){
        var segs_end=document.cookie.indexOf(';',segs_beg);
        if(segs_end==-1)segs_end=document.cookie.length;
        rsi_segs=document.cookie.substring(segs_beg,segs_end).split('|');
    }
}
var segLen=20
var segQS="",segArr=new Array()
if (rsi_segs.length<segLen) {
    segLen=rsi_segs.length
}
for (var i=0;i<segLen;i++){
    segArr=rsi_segs[i].split("_")
    if (segArr.length>1) segQS+=("rsi"+"="+segArr[1]+"+") 
}