'use strict';

var Parse = require('parse/node');

var parse_service = require('./parse_service');
// var crud_service = require('./crud_service');
// var race_crawler = require('./race_crawler');
var rest_request_utils = require("../utils/rest_request_utils");
var string_utils = require("../utils/string_utils");
var racecard_parser = require("../parsers/racecard_parser");
var url_util = require('url');
var _ = require('underscore');
var async = require('async');

// var RacecardModel = require("../models/racecard");

var async_callback = function(err, results){
    if(err){
        console.log("Error " +err);
    } else {
        console.log("Success   " + result);
    }
};

module.exports = {
    racecard: function(url, options){
        var url_parts = url_util.parse(url, true);
        rest_request_utils.get(url,null,function(html){

            var htmlMD5 = string_utils.hash_string(html);

            ////////////////***************
            var racecard_json = racecard_parser.parse(html, {
                dt:url_parts.path.split('/')[3],
                page_url:url,
                page_hash:htmlMD5
            });

            _.forEach(racecard_json.meetings,function(meeting, meeting_index, meeting_list){
                _.forEach(meeting.races,function(race, race_index, race_list){
                    meeting.races[race_index] = parse_service.to_parse_obj('Race',race);
                });
                racecard_json.meetings[meeting_index] = parse_service.to_parse_obj('Meeting', meeting);
            });
            racecard_json = parse_service.to_parse_obj('Racecard',racecard_json);

            var tmp_list = []
            _.forEach(racecard_json.get('meetings'),function(meeting, meeting_index, meeting_list){
                tmp_list.push(meeting.get('races'));
            });

            Parse.Object.saveAll(tmp_list)
            .then(
                function(obj){
                    Parse.Object.saveAll(racecard_json.get('meetings'))
                        .then(
                            function(obj1){
                                racecard_json.save(null, {
                                    success: function(result) {
                                        // Execute any logic that should take place after the object is saved.
                                        console.log('New object created with objectId: ' + result.id);
                                    },
                                    error: function(result, error) {
                                        // Execute any logic that should take place if the save fails.
                                        // error is a Parse.Error with an error code and message.
                                        console.log('Failed to create new object, with error code: ' + error.message);
                                    }
                                });
                            });
                    return 'saved';
                }
            );

            ////////////////***************
            
            // RacecardModel.findByPageUrl(url, function(error, racecards){
            //     if(!error){
            //         var racecard = racecards[0];
            //         if (racecard !== undefined && racecard.page_body_hash === htmlMD5){
            //             console.log(url + " body has not changed, doing nothing");
            //         } else {
            //
            //             //parse doc
            //             var racecard_json = racecard_parser.parse(html, {
            //                 dt:url_parts.path.split('/')[3],
            //                 page_url:url,
            //                 page_hash:htmlMD5
            //             });
            //
            //             //set new docment
            //             var racecardNew = new RacecardModel(racecard_json);
            //
            //             //save new document
            //             crud_service.save(racecardNew);
            //             //remove previous
            //             if(racecard !== undefined) {
            //                 racecard.remove(function (error) {
            //                     if (error) {
            //                         console.log("An error happened -> " + JSON.stringify(error));
            //                     }
            //                     console.log("Document was removed!");
            //                 });
            //             }
            //
            //             racecardNew.cards.forEach(function(meeting, index, array){
            //                 race_crawler.race(meeting);
            //             });
            //             console.log(url + " body changed, updating");
            //         }
            //     } else {
            //         console.log(JSON.stringify(error));
            //     }
            // });
            //

        });
    },

}