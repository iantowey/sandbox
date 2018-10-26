'use strict';

var couchbase_service = require("../services/couchbase_service");
var ottoman = require("ottoman");

ottoman.bucket = couchbase_service.get_bucket();

var RacecardMdl = ottoman.model("Racecard", {
    createdON: { type: "Date", default: function() { return new Date() } },
    date:"string",
    page_url: "string",
    page_body_hash: "string",
    cards: [{
        track: "string",
        track_url:"string",
        track_going:"string",
        track_surface:"string",
        races:[{
            race_time:"string",
            race_name:"string",
            racecard_url:"string",
            result_url:"string",
            abandoned:"string"
        }]
    }]
}, {
    index: {
        findByPageUrl: {
            by: "page_url",
            type: "n1ql"
        },
        findByDate: {
            by: "date",
            type: "n1ql"
        }
    }
}
);

ottoman.ensureIndices(function(err) {
    if (err) return console.error(err);
});

module.exports=RacecardMdl;