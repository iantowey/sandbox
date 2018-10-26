'use strict';

var cheerio = require('cheerio');
var crypto = require('crypto');
var string_utils = require('./../../utils/string_utils')
var md5sum = crypto.createHash('md5');

module.exports = {
    parse: function(html,options){
        var $ = cheerio.load(html);
        md5sum.update($('html').toString());

        var daily_racecard = {
            date:options.dt,
            crawl_timestamp:new Date().getTime(),
            page_url:options.page_url,
            page_body_hash:md5sum.digest('hex'),
            cards:[]
        };

        $('.r .rac-dgp').each(function(i, element){
            var race_card_section = $(this);
            var racecard = {
                track:race_card_section.find('.hdr.t2 a').text().trim(),
                track_url:race_card_section.find('.hdr.t2 a').attr('href'),
                track_going:race_card_section.find('.list.hdr-sublinks li:nth-child(1):has(strong)').html().replace(/<strong>.*<\/strong>/g,''),
                track_surface:string_utils.nvl(race_card_section.find('.list.hdr-sublinks li:nth-child(2):has(strong)').html()).replace(/<strong>.*<\/strong>/g,''),
                races:[]
            };
            var tmp_races = [];
            race_card_section.find('.rac-cards').each(function(i, element){
                var race_section = $(this);
                tmp_races.push({
                    race_time:race_section.find('.ix.ixt').text(),
                    race_name:race_section.find('.ix.ixc').text().replace(/\n/g,' ').replace(/\s+/g,' '),
                    racecard_url:race_section.find('.ix.ixc a').attr('href'),
                    result_url:race_section.find('.ix.ixv a').attr('href')
                })
                //time
            });
            racecard.races = tmp_races;
            daily_racecard.cards.push(racecard);
        });

        daily_racecard = string_utils.scrub_obj_whitespace(daily_racecard);

        return daily_racecard;
    }
}


