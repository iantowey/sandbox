'use strict';

var cheerio = require('cheerio');
var string_utils = require('../utils/string_utils');

module.exports = {
  parse: function(html, options){

      
      var race = {
          page_url:options.page_url,
          page_body_hash:options.page_body_hash,
          race_name:null,
          race_date:options.race_date,
          race_time:null,
          race_distance:null,
          race_class:null,
          prize_money:null,
          race_going:null,
          race_surface:null,
          runners:null
      }
      
      var race_runners = []
      var $ = cheerio.load(html);


      $('#racecard tbody tr:not(.disabled)').each(function(i, element){
          var race_entry = $(this);
          race_runners.push({
              running: true,
              stall:race_entry.find('td:nth-child(1)').text().replace(/(\n|\s+)/g,' ').trim(),
              form:race_entry.find('td:nth-child(2)').text().replace(/(\n|\s+)/g,' ').trim(),
              silk_img_link:race_entry.find('td:nth-child(3) img').attr('src'),
              horse:race_entry.find('td:nth-child(3) a').text().replace(/(\n|\s+)/g,' ').trim(),
              profile_url:race_entry.find('td:nth-child(3) a').attr('href'),
              breeding:race_entry.find('td:nth-child(3) a').attr('title'),
              age:race_entry.find('td:nth-child(4)').text(),
              weight:race_entry.find('td:nth-child(5)').text(),
              trainer:{
                  name:race_entry.find('td:nth-child(6) a').text(),
                  profile_url:race_entry.find('td:nth-child(6) a').attr('href')
              },
              jockey:{
                  name:race_entry.find('td:nth-child(7) a').text(),
                  profile_url:race_entry.find('td:nth-child(7) a').attr('href')
              },
              or:race_entry.find('td:nth-child(8)').text(),
              sp:race_entry.find('td:nth-child(10)').text()
          });
      });

      //non runners
      $('#racecard tbody tr.disabled').each(function(i, element){
          var race_entry = $(this);
          race_runners.push({
              running: false,
              stall:race_entry.find('td:nth-child(1)').text().replace(/(\n|\s+)/g,' ').trim(),
              form:race_entry.find('td:nth-child(2)').text().replace(/(\n|\s+)/g,' ').trim(),
              silk_img_link:race_entry.find('td:nth-child(3) img').attr('src'),
              horse:race_entry.find('td:nth-child(3) strong').text().replace(/(\n|\s+)/g,' ').trim(),
              profile_url:race_entry.find('td:nth-child(3) a').attr('href'),
              breeding:race_entry.find('td:nth-child(3) a').attr('title'),
              age:race_entry.find('td:nth-child(4)').text(),
              weight:race_entry.find('td:nth-child(5)').text(),
              trainer:{
                  name:race_entry.find('td:nth-child(6) a').text(),
                  profile_url:race_entry.find('td:nth-child(6) a').attr('href')
              },
              jockey:{
                  name:race_entry.find('td:nth-child(7) a').text(),
                  profile_url:race_entry.find('td:nth-child(7) a').attr('href')
              },
              or:race_entry.find('td:nth-child(8)').text(),
              sp:race_entry.find('td:nth-child(10)').text()
          });
      });
      race.runners = string_utils.scrub_obj_whitespace(race_runners);

      return race;
  }
};

