'use strict';

var string_utils = require('./../../utils/string_utils');
var cheerio = require('cheerio');

module.exports = {
    parse:function(html, options){
        var horse_details = {};

        var $ = cheerio.load(html);

        var racecard_header = $('.racecard-header.t3 .content-header');

        horse_details = {
            name : racecard_header.find('h2').text(),
            age : racecard_header.find('ul li:nth-child(1)').text(),
            sex :  racecard_header.find('ul li:nth-child(2)').text(),
            profile_url: (('url' in options) ? options['url'] : null),
            sire : {
                name:racecard_header.find('ul li:nth-child(3) a').text(),
                profile_url:racecard_header.find('ul li:nth-child(3) a').attr('href')
            },
            dam :  {
                name:racecard_header.find('ul li:nth-child(4) a').text(),
                profile_url:racecard_header.find('ul li:nth-child(4) a').attr('href')
            },
            trainer :  {
                name:racecard_header.find('ul li:nth-child(6) a:nth-child(2)').text(),
                profile_url:racecard_header.find('ul li:nth-child(6) a:nth-child(2)').attr('href')
            },
            owner :  racecard_header.find('ul li:nth-child(7)').text(),
            historical_form : []
        };

        var tmp_hist_summary = $('.tbl.t3:not(#horse-profile-form) tbody tr');
        horse_details.historical_summary = {
            period:tmp_hist_summary.find('td:nth-child(1)').text(),
            runs:tmp_hist_summary.find('td:nth-child(2)').text(),
            first:tmp_hist_summary.find('td:nth-child(3)').text(),
            second:tmp_hist_summary.find('td:nth-child(4)').text(),
            third:tmp_hist_summary.find('td:nth-child(5)').text(),
            win_percentage:tmp_hist_summary.find('td:nth-child(6)').text(),
            level_stake:tmp_hist_summary.find('td:nth-child(7)').text()
        };

        $('#horse-profile-form tbody tr').each(function(i, element){
            var lifetime_stat = $(this);
            horse_details.historical_form.push({
                dt:lifetime_stat.find('td:nth-child(1) a').text(),
                race_notes:lifetime_stat.find('td:nth-child(1) a').attr('title'),
                finish_position:lifetime_stat.find('td:nth-child(2)').text(),
                weight:lifetime_stat.find('td:nth-child(3)').text(),
                bha:lifetime_stat.find('td:nth-child(4)').text(),
                race_details:lifetime_stat.find('td:nth-child(5)').text(),
                jockey:{
                    name:lifetime_stat.find('td:nth-child(6) a').text(),
                    profile_url:lifetime_stat.find('td:nth-child(6) a').attr('href')
                },
                trainer:{
                    name:lifetime_stat.find('td:nth-child(7) a').text(),
                    profile_url:lifetime_stat.find('td:nth-child(7) a').attr('href')
                },
                distance_back:lifetime_stat.find('td:nth-child(9)').text(),
                sp:lifetime_stat.find('td:nth-child(10)').text()
            }) ;
        });

        horse_details = string_utils.scrub_obj_whitespace(horse_details);

        return horse_details;
    }
};