var chai   = require("chai");
var horse_details_parser = require("../services/crawlers/parsers/horse_details_parser");
var fs = require('fs');

var file_name = './tmp/apache-outlaw.html';

describe("Parse horse details", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });
    
    it("get horse details", function (done) {

        var html = fs.readFileSync(file_name);
        var horse_details = horse_details_parser.parse(html);

        console.log(JSON.stringify(horse_details));
        chai.expect(horse_details.name).equal('Apache Outlaw (IRE)');
        chai.expect(horse_details.sex).equal('Sex: Bay Gelding');
        chai.expect(horse_details.sire.name.trim()).equal('Westerner');
        chai.expect(horse_details.dam.name.trim()).equal('Bermuda Bay (IRE)');
        chai.expect(horse_details.trainer.name.trim()).equal('Miss R Curtis');
        chai.expect(horse_details.owner).equal('Owner: Mr G Costelloe');
        chai.expect(horse_details.historical_form.length).to.equal(5);
        done();
    });
});