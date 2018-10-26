'use strict';

var couchbase_service = require("../services/couchbase_service");
var ottoman = require("ottoman");

ottoman.bucket = couchbase_service.get_bucket();

var ResultMdl = ottoman.model("Result", {
        createdON: { type: "Date", default: function() { return new Date() } },
        page_url: "string",
        page_body_hash: "string",
        race_name:"string",
        race_date:"string",
        race_time:"string",
        race_distance:"string",
        race_class:"string",
        classification:"string",
        prize_money:"string",
        race_going:"string",
        race_surface:"string",
        winning_time:"string",
        start_time:"string",

        results:[{
            position: "string",
            distance: "string",
            weight: "string",
            horse:{
                no:"string",
                name: "string",
                profile_url: "string",
                age: "string"
            },
            trainer: {
                name: "string",
                profile_url: "string"
            },
            jockey: {
                name: "string",
                profile_url: "string"
            },
            or: "string",
            sp: "string"
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

module.exports=ResultMdl;


