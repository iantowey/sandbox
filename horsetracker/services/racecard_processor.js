'use strict';

var fs = require('fs');
var cheerio = require('cheerio');
var dateFormat = require('dateformat');
var crypto = require('crypto');

var kafka = require('kafka-node'),
    Producer = kafka.Producer,
    KeyedMessage = kafka.KeyedMessage,
    client = new kafka.Client(),
    producer = new Producer(client);

var md5sum = crypto.createHash('md5');

var today_dt = dateFormat(new Date(), "dd-mm-yyyy");
var file_name = '../tmp/racecard_' + today_dt + ".html";

var $ = cheerio.load(fs.readFileSync(file_name));
md5sum.update($('html').toString());

var daily_racecard = {
    date:today_dt,
    crawl_timestamp:new Date().getTime(),
    page_url:'http://www.sportinglife.com/racing/racecards/' + today_dt,
    page_body_hash:md5sum.digest('hex'),
    cards:[]
};

$('.r .rac-dgp').each(function(i, element){
    var race_card_section = $(this);
    var racecard = {
        track:race_card_section.find('.hdr.t2 a').text().trim(),
        track_url:race_card_section.find('.hdr.t2 a').attr('href'),
        track_going:race_card_section.find('.list.hdr-sublinks li:nth-child(1)').text(),
        track_surface:race_card_section.find('.list.hdr-sublinks li:nth-child(2)').text(),
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

console.log(JSON.stringify(daily_racecard));

