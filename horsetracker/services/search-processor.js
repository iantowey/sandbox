'use strict';

var fs = require('fs');
var cheerio = require('cheerio');

var file_name = '../tmp/search-result.html';

var $ = cheerio.load(fs.readFileSync(file_name));

var search_result = {
    query_string:'star',
    type:'H',
    crawl_timestamp:new Date().getTime(),
    page_url:'http://www.sportinglife.com/ajax/racing/profile-search?search=star&type=H&page=1',
    results:[]
};

var num_result_str = $('.content-header.alt p').text();
var total_result_count = num_result_str.match(/[0-9]* (?=results)/g)[0].trim();

console.log(num_result_str);
console.log(total_result_count);

$('.search-results .search-item.t2').each(function(i, element){     
    var search_result_tmp = $(this);
    search_result.results.push({
       name: search_result_tmp.find('h4').text().trim().replace(/\s+/g,' '),
        profile_url: search_result_tmp.find('a').attr('href')
    });
});

console.log(JSON.stringify(search_result));

//$('.r .rac-dgp').each(function(i, element){
//    var race_card_section = $(this);
//    var racecard = {
//        track:race_card_section.find('.hdr.t2 a').text().trim(),
//        track_url:race_card_section.find('.hdr.t2 a').attr('href'),
//        track_going:race_card_section.find('.list.hdr-sublinks li:nth-child(1)').text(),
//        track_surface:race_card_section.find('.list.hdr-sublinks li:nth-child(2)').text(),
//        races:[]
//    };
//    var tmp_races = [];
//    race_card_section.find('.rac-cards.fr-complete.click').each(function(i, element){
//        var race_section = $(this);
//        //time race_section.find('.ix.ixt').text()
//        console.log('here');
//    });
//    racecard.races = tmp_races;
//    daily_racecard.cards.push(racecard);
//});

//console.log(JSON.stringify(daily_racecard));

