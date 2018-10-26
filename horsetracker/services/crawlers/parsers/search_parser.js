'use strict';

var cheerio = require('cheerio');
var html_utils = require("../../utils/html_utils");

module.exports = {
    parse: function(html, options){
        var $ = cheerio.load(html);

        var search_result = {
            crawl_timestamp:options.crawl_timestamp,
            url_details:html_utils.parse_url(options.page_url),
            result_count:null,
            results:[]
        };

        var num_result_str = $('.content-header.alt p').text();
        search_result.result_count = num_result_str.match(/[0-9]* (?=results)/g)[0].trim();

        $('.search-results .search-item.t2').each(function(i, element){
            var search_result_tmp = $(this);
            search_result.results.push({
                name: search_result_tmp.find('h4').text().trim().replace(/\s+/g,' '),
                profile_url: search_result_tmp.find('a').attr('href')
            });
        });

        return search_result;
    }
}

