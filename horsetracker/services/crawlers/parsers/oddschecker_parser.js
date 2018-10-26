'use strict';

var html_utils = require("./../..//utils/html_utils");
var string_utils = require('./../../utils/string_utils');
var cheerio = require('cheerio');

module.exports = {
    parse: function (html, options) {
        var oddschecker_data = {
            race_info:html_utils.parse_url(options.url),
            race_details:null,
            bookmakers: [],
            entries_odds:[]
        };

        var $ = cheerio.load(html);

        oddschecker_data.race_details = $('.page-description.module .event p').text();

        var oddschecker_tbl = $('#oddsTableContainer table');

        //find bookmakers
        oddschecker_tbl.find('thead .eventTableHeader td[data-bk]').each(function (i, element) {
            var header_row_col = $(this);
            oddschecker_data.bookmakers.push(header_row_col.find('a').attr('title'));
        });

        //find odds
        oddschecker_tbl.find('tbody tr').each(function (i, element) {
            var body_row = $(this);
            var entry_odds = {
                stall:body_row.find('td:nth-child(1)').text(),
                form:body_row.find('td:nth-child(2)').text(),
                horse:{
                    name:body_row.find('td:nth-child(3) .selTxt').text(),
                    profile:body_row.find('td:nth-child(2) a').attr('href'),
                    race_odds_history:body_row.find('td:nth-child(3) a:nth-child(1)').attr('href'),
                    silk_url:body_row.find('td:nth-child(3) img').attr('src')
                },
                odds:[]
            }

            body_row.find('td[data-hcap]').each(function(i, element){
                var row_col = $(this);
                entry_odds.odds.push({
                    fractional:row_col.text(),
                    decimal:row_col.attr('data-odig')
                })
            });
            oddschecker_data.entries_odds.push(entry_odds);
        });

        oddschecker_data = string_utils.scrub_obj_whitespace(oddschecker_data);

        return oddschecker_data;
    }
}