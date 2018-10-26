var chai   = require("chai");
var race_entries_parser = require("../services/crawlers/parsers/race_entries_parser");
var fs = require('fs');

var file_name = './tmp/execter-30-03-2016.html';

describe("Parse horse details", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });

    it("get horse details", function (done) {

        var html = fs.readFileSync(file_name);
        var race_entries = race_entries_parser.parse(html);

        console.log(JSON.stringify(race_entries));
        chai.expect(race_entries[0].running).to.equal(true);
        chai.expect(race_entries[0].stall).to.equal('1');
        chai.expect(race_entries[0].form).to.equal('F-2517');
        chai.expect(race_entries[0].silk_img_link).to.equal('http://images.sportinglife.com/racing/tsilks/127652.png');
        chai.expect(race_entries[0].horse).to.equal('Apache Outlaw (IRE)');
        chai.expect(race_entries[0].profile_url).to.equal('/racing/profiles/horse/914744/apache-outlaw');
        done();
    });
});