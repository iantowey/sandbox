var html_utils = require("../services/utils/html_utils");

describe("Parse url", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });

    it("parse url", function (done) {

        var url = "http://www.sportinglife.com/ajax/racing/profile-search?search=star&type=H&page=1";
        var url_parts = html_utils.parse_url(url);

        console.log(JSON.stringify(url_parts));
        done();
    });
});