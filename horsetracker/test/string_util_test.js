var string_utils = require("../services/utils/string_utils");
var fs = require('fs');



describe("String util test", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });

    it("Bug fix for : Cannot assign to read only property '0' of Bet 365", function (done) {

        var file_name = './tmp/oddschecker_json_err';
        var json = JSON.parse(fs.readFileSync(file_name));

        console.log(JSON.stringify(string_utils.scrub_obj_whitespace(json)));
        done();
    });

    it("Bug fix for : clean string", function (done) {
        var file_name = './tmp/detail_json_test';
        var json = JSON.parse(fs.readFileSync(file_name));

        console.log(JSON.stringify(string_utils.scrub_obj_whitespace(json)));
        done();
    });
});