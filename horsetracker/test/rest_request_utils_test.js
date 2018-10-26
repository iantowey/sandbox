var chai   = require("chai");
var rest_request_utils = require("../services/utils/rest_request_utils");
var date_utils = require("../services/utils/date_utils");


describe("GET url request tests", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });

    describe("GET todays racecard", function() {
        it("GET todays racecard 200", function (done) {
            var url = 'http://www.sportinglife.com/racing/racecards/' + date_utils.today_yy_mm_yyyy();
            rest_request_utils.get(url, null, function (data) {
                chai.expect(data).to.contain('</html>');
                done();
            });
        });
    });
});