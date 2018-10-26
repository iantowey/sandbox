var chai   = require("chai");
var oddschecker_parser = require("../services/crawlers/parsers/oddschecker_parser");
var fs = require('fs');

var file_name = './tmp/Oddschecker.html';

describe("Parse oddschecker race", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });
    
    it("Parse oddschecker race tesr", function (done) {

        var html = fs.readFileSync(file_name);
        var oddschecker_json = oddschecker_parser.parse(html,{url:"http://www.oddschecker.com/horse-racing/2016-04-05-newton-abbot/13:50/winner"});

        console.log(JSON.stringify(oddschecker_json));
        // chai.expect(horse_details.name).equal('Apache Outlaw (IRE)');
        // chai.expect(horse_details.sex).equal('Sex: Bay Gelding');
        // chai.expect(horse_details.sire.name).equal('Westerner');
        // chai.expect(horse_details.dam.name).equal('Bermuda Bay (IRE)');
        // chai.expect(horse_details.trainer.name).equal('Miss R Curtis');
        // chai.expect(horse_details.owner).equal('Owner: Mr G Costelloe');
        // chai.expect(horse_details.historical_form.length).to.equal(5);
        done();
    });
});