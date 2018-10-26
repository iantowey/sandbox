'use strict';

var cheerio = require('cheerio');
var string_utils = require('./../../utils/string_utils');

module.exports = {

  parse: function(html, options){

      var race_result = {}
      var $ = cheerio.load(html);

      var racecard_header = $('.racecard-header .content-header');

      race_result.race_name = racecard_header.find('h2').html().replace(/<span.*<\/span>/g,'') ;
      race_result.race_time = $('#content .header-nav h1').text() ;
      race_result.classification = racecard_header.find('ul li:nth-child(1)').text();
      race_result.prize_money = racecard_header.find('ul li:nth-child(2)').text();
      race_result.going = racecard_header.find('ul li:nth-child(3)').text();
      race_result.surface = racecard_header.find('ul li:nth-child(4)').text();
      race_result.winning_time = $('.racecard-header .racecard-status ul li:nth-child(1)').text();
      race_result.start_time = $('.racecard-header .racecard-status ul li:nth-child(2)').text();

      race_result.results = [];
      $('#racecard-tabs .tbl.t3 tbody tr').each(function(i, element){
          var row = $(this);
          race_result.results.push({
              position:row.find('td:nth-child(1)').text(),
              distance:row.find('td:nth-child(2)').text(),
              name:row.find('td:nth-child(3) a').text(),
              horse_profile_url:row.find('td:nth-child(3) a').attr('href'),
              trainer:row.find('td:nth-child(4) a').text(),
              trainer_profile_url:row.find('td:nth-child(4) a').attr('href'),
              age:row.find('td:nth-child(5)').text(),
              weight:row.find('td:nth-child(6)').text(),
              jockey:row.find('td:nth-child(7) a').text(),
              jockey_profile_url:row.find('td:nth-child(7) a').attr('href'),
              sp:row.find('td:nth-child(8)').text(),

          });
      });

      race_result = string_utils.scrub_obj_whitespace(race_result)
      return race_result;
  }
};

