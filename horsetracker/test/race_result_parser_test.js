var chai   = require("chai");
var race_result_parser = require("../services/crawlers/parsers/race_result_parser");
var fs = require('fs');

var file_name = './tmp/result.html';

describe("Parse horse details", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });

    it("get horse details", function (done) {

        var html = fs.readFileSync(file_name);
        var race_result = race_result_parser.parse(html);

        console.log(JSON.stringify(race_result));
        chai.expect(race_result.race_name.trim()).to.equal('Betway Spring Mile (Handicap) (Str)');
        chai.expect(race_result.race_time.trim()).to.equal('14:15 Doncaster');
        chai.expect(race_result.going.trim()).to.equal('Going: Soft');
        chai.expect(race_result.surface.trim()).to.equal('Surface: Turf');
        chai.expect(race_result.results.length).to.equal(22);

        done();
    });
});