'use strict';

var cheerio = require('cheerio');
var string_utils = require('../utils/string_utils');

module.exports = {

    parse: function(html, options){

        var race_result = {}
        var $ = cheerio.load(html);

        var racecard_header = $('.racecard-header .content-header');

        race_result.page_url = options.page_url;
        race_result.page_body_hash = options.page_body_hash;
        race_result.race_name = racecard_header.find('h2').html().replace(/<span.*<\/span>/g,'') ;
        race_result.race_time = $('#content .header-nav h1').text() ;
        race_result.classification = racecard_header.find('ul li:nth-child(1)').text();
        race_result.prize_money = racecard_header.find('ul li:nth-child(2)').text();
        race_result.race_going = racecard_header.find('ul li:nth-child(3)').text();
        race_result.race_surface = racecard_header.find('ul li:nth-child(4)').text();
        race_result.winning_time = $('.racecard-header .racecard-status ul li:nth-child(1)').text();
        race_result.start_time = $('.racecard-header .racecard-status ul li:nth-child(2)').text();

        race_result.results = [];
        $('#racecard-tabs .tbl.t3 tbody tr:not(.note)').each(function(i, element){
            var row = $(this);
            race_result.results.push({
                position:row.find('td:nth-child(1)').text(),
                distance:row.find('td:nth-child(2)').text(),
                weight:row.find('td:nth-child(6)').text(),
                horse: {
                    no: row.find('td:nth-child(3) .name').html().replace(/<a(.|\n)*/g, "").trim(),
                    name: row.find('td:nth-child(3) a').text(),
                    profile_url: row.find('td:nth-child(3) a').attr('href'),
                    age:row.find('td:nth-child(5)').text(),
                },
                trainer:{
                    name: row.find('td:nth-child(4) a').text(),
                    profile_url: row.find('td:nth-child(4) a').attr('href'),
                },
                jockey: {
                    name: row.find('td:nth-child(7) a').text(),
                    profile_url: row.find('td:nth-child(7) a').attr('href'),
                },
                sp:row.find('td:nth-child(8)').text(),

            });
        });

        race_result = string_utils.scrub_obj_whitespace(race_result)
        return race_result;
    }
};
