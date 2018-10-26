'use strict';

var fs = require('fs');
var cheerio = require('cheerio');

var file_name = '../tmp/execter-30-03-2016.html';

var $ = cheerio.load(fs.readFileSync(file_name));

var race_entries = []
$('#racecard tbody tr:not(.disabled)').each(function(i, element){
    var race_entry = $(this);
    race_entries.push({
        running: true,
        stall:race_entry.find('td:nth-child(1)').text().replace(/(\n|\s+)/g,' ').trim(),
        form:race_entry.find('td:nth-child(2)').text().replace(/(\n|\s+)/g,' ').trim(),
        silk_img_link:race_entry.find('td:nth-child(3) img').attr('src'),
        horse:race_entry.find('td:nth-child(3) a').text().replace(/(\n|\s+)/g,' ').trim(),
        horse_link:race_entry.find('td:nth-child(3) a').attr('href')
        // ,
        //
        // a:race_entry.find('td:nth-child(2)').text(),
        // weight:race_entry.find('td:nth-child(2)').text(),
        // trainer:race_entry.find('td:nth-child(2)').text(),
        // jockey:race_entry.find('td:nth-child(2)').text(),
        // or:race_entry.find('td:nth-child(2)').text(),
        // naps:race_entry.find('td:nth-child(2)').text(),
        // sp:race_entry.find('td:nth-child(2)').text()
    });
});
//non runners
$('#racecard tbody tr.disabled').each(function(i, element){
    var race_entry = $(this);
    race_entries.push({
        running: false,
        stall:race_entry.find('td:nth-child(1)').text().replace(/(\n|\s+)/g,' ').trim(),
        form:race_entry.find('td:nth-child(2)').text().replace(/(\n|\s+)/g,' ').trim(),
        silk_img_link:race_entry.find('td:nth-child(3) img').attr('src'),
        horse:race_entry.find('td:nth-child(3) strong').text().replace(/(\n|\s+)/g,' ').trim(),
        horse_link:race_entry.find('td:nth-child(3) a').attr('href')
        // ,
        //
        // a:race_entry.find('td:nth-child(2)').text(),
        // weight:race_entry.find('td:nth-child(2)').text(),
        // trainer:race_entry.find('td:nth-child(2)').text(),
        // jockey:race_entry.find('td:nth-child(2)').text(),
        // or:race_entry.find('td:nth-child(2)').text(),
        // naps:race_entry.find('td:nth-child(2)').text(),
        // sp:race_entry.find('td:nth-child(2)').text()
    });
});

console.log(JSON.stringify(race_entries));