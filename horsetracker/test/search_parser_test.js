var fs = require('fs');
var chai   = require("chai");
var search_parser = require("../services/crawlers/parsers/search_parser");
var file_name = './tmp/search-result.html';


describe("Parse search results", function() {

    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });
    
    it("get search results", function (done) {

        var html = fs.readFileSync(file_name);
        var search_result = search_parser.parse(html,{
            page_url:'http://www.sportinglife.com/ajax/racing/profile-search?search=star&type=H&page=1',
            crawl_timestamp:new Date().getTime()
        });

        console.log(JSON.stringify(search_result));
        chai.expect(search_result.url_details.query.search).equal('star');
        chai.expect(search_result.url_details.query.type).equal('H');
        chai.expect(search_result.url_details.host).equal('www.sportinglife.com');
        chai.expect(search_result.result_count).equal('8764');
        chai.expect(search_result.results.length).equal(10);
        done();
    });
});