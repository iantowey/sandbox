var RaceModel = require("../models/race");
var _ = require('underscore');
var neo4j_service = require('../services/neo4j_service');
var dateformat = require('dateformat');
var async = require('async');

String.prototype.clean_str = function(){
    return this.trim().replace(/ /g,"_").replace(/,/g,"");
}

var today_yy_mm_yyyy =  function() {
    return dateformat(new Date(), "dd-mm-yyyy");
}

var neo_callback = function(err, results) {
    if (err) {
        console.log("ERROR: create response  " + JSON.stringify(err) );
    }
    console.log("Success response  " + JSON.stringify(results));
}

module.exports = {

    build_graph : function (race_url) {
        
        console.log('Processing ' + race_url);
        RaceModel.findByPageUrl(race_url, function (error, result) {
            if (error) {
                console.log("An error happened -> " + JSON.stringify(error));
            }
            var race = result[0];
            _.each(race.runners, function (element, index, list) {

                var entry = element;

                neo4j_service.create_horse({
                    name: entry.horse.clean_str(),
                    profile_url: entry.profile_url.trim(),
                    dt: today_yy_mm_yyyy()
                }, neo_callback);

                neo4j_service.create_trainer({
                    name: entry.trainer.name.trim().clean_str(),
                    profile_url: entry.trainer.profile_url,
                    dt: today_yy_mm_yyyy()
                }, neo_callback);

                neo4j_service.create_jockey({
                    name: entry.jockey.name.clean_str(),
                    profile_url: entry.jockey.profile_url,
                    dt: today_yy_mm_yyyy()
                }, neo_callback);
            });

            setTimeout(function () {
                _.each(race.runners, function (element, index, list) {
                    var entry = element;

                    async.series([

                        function (callback) {
                            neo4j_service.create_relationship_horse_trainer({
                                trainer_name: entry.trainer.name.clean_str(),
                                runner_name: entry.horse.clean_str(),
                                rel_dt: today_yy_mm_yyyy()
                            }, neo_callback);
                            callback();
                        },
                        function (callback) {
                            neo4j_service.create_relationship_horse_jockey({
                                jockey_name: entry.jockey.name.clean_str(),
                                runner_name: entry.horse.clean_str(),
                                rel_dt: today_yy_mm_yyyy()
                            }, neo_callback);
                            callback();
                        }
                    ], function (err, results) {
                        if (err)
                            console.log(JSON.stringify(err));
                        console.log(JSON.stringify(results));
                    });
                });
                console.log("All Done");
            }, 6000);
        });
    }
}