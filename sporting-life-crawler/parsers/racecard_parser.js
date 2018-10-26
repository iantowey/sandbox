'use strict';

var parse_service = require('../services/parse_service');
var cheerio = require('cheerio');
var string_utils = require('../utils/string_utils')

module.exports = {
    parse: function(html,options){
        var $ = cheerio.load(html);

        var daily_racecard = {
            createdON:new Date(),
            date:options.dt,
            // crawl_timestamp:new Date().getTime(),
            page_url:options.page_url,
            page_body_hash:options.page_hash,
            meetings:[]
        };

        $('.r .rac-dgp').each(function(i, element){
            var race_card_section = $(this);
            var meeting = {
                track:race_card_section.find('.hdr.t2 a').text().trim(),
                track_url:race_card_section.find('.hdr.t2 a').attr('href'),
                track_going:race_card_section.find('.list.hdr-sublinks li:nth-child(1):has(strong)').html().replace(/<strong>.*<\/strong>/g,''),
                track_surface:string_utils.nvl(race_card_section.find('.list.hdr-sublinks li:nth-child(2):has(strong)').html()).replace(/<strong>.*<\/strong>/g,''),
                races:[]
            };
            var tmp_races = [];
            race_card_section.find('.rac-cards').each(function(i, element){
                var race_section = $(this);
                var tmp_race = {
                    race_time:race_section.find('.ix.ixt').text(),
                    race_name:race_section.find('.ix.ixc').text().replace(/\n/g,' ').replace(/\s+/g,' '),
                    racecard_url:race_section.find('.ix.ixc a').attr('href'),
                    result_url:race_section.find('.ix.ixv a').attr('href'),
                    abandoned:race_section.find('.ix.ixv .dnf').text()
                };
                //var race_parse = parse_service.save_obj('Race',tmp_race);
                tmp_races.push(tmp_race);
                //time
            });
            meeting.races = tmp_races;
            //var meeting_parse = parse_service.save_obj('Meeting',meeting);
            daily_racecard.meetings.push(meeting);
        });

        // daily_racecard = string_utils.scrub_obj_whitespace(daily_racecard);

        return daily_racecard;
    }
}


