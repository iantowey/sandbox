var chai   = require("chai");
var racecard_parser = require("../services/crawlers/parsers/racecard_parser");
var fs = require('fs');
var dateFormat = require('dateformat');

var today_dt = dateFormat(new Date(), "dd-mm-yyyy");
// var file_name = './tmp/racecard_' + today_dt + '.html';
var file_name = './tmp/racecard_01-04-2016.html';

describe("Parse racecard details", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });

    it("get todays racecards details", function (done) {

        var html = fs.readFileSync(file_name);
        var racecards = racecard_parser.parse(html,{dt:today_dt, page_url:'http://www.sportinglife.com/racing/racecards/' + today_dt});

        console.log(JSON.stringify(racecards));
        chai.expect(racecards.cards.length).to.equal(16);
        chai.expect(racecards.cards[0].track).to.equal('Newbury');
        chai.expect(racecards.cards[0].track_surface.trim()).to.equal('Turf');
        chai.expect(racecards.cards[0].races.length).to.equal(7);
        chai.expect(racecards.cards[0].races[0].race_time).to.equal("13:50");
        chai.expect(racecards.cards[0].races[0].race_name).to.contain(" Premier Novices' Hurdle (4yo+, Class 3, 2m 4f 118y,");
        done();
    });
});