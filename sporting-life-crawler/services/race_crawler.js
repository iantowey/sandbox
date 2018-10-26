/**
 * Created by ian on 10/04/16.
 */

var conf = require('../config/' + process.env.NODE_ENV + '.json');

var crud_service = require('./crud_service');
var rest_request_utils = require("../utils/rest_request_utils");
var race_parser = require('../parsers/race_parser');
var result_parser = require('../parsers/result_parser');
var string_utils = require('../utils/string_utils');
var url_util = require('url');

var _ = require('underscore');
var RaceModel = require('../models/race');
var ResultModel = require('../models/result');

var get_runners = function(url){

    var url_parts = string_utils.url_parse(url);
    rest_request_utils.get(url,null,function(html){

        var htmlMD5 = string_utils.hash_string(html);
        RaceModel.findByPageUrl(url, function(error, races){
            if(!error){
                var race = races[0];
                if (race !== undefined && race.page_body_hash === htmlMD5){
                    console.log(url + " already crawled, doing nothing");
                } else {
                    //parse doc
                    var race_json = race_parser.parse(html, {
                        page_url:url,
                        page_body_hash:htmlMD5,
                        race_date:url_parts.path_array[3]
                    });

                    //set new docment
                    var raceNew = new RaceModel(race_json);

                    //save new document
                    crud_service.save(raceNew);
                }
            } else {
                console.log(JSON.stringify(error));
            }
        });
    });
}

var get_result = function(url){
    var url_parts = url_util.parse(url, true);

    rest_request_utils.get(url,null,function(html){

        var htmlMD5 = string_utils.hash_string(html);
        ResultModel.findByPageUrl(url, function(error, results){
            if(!error){
                var result = results[0];
                if (result !== undefined && result.page_body_hash === htmlMD5){
                    console.log(url + " already crawled, doing nothing");
                } else {
                    //parse doc
                    var result_json = result_parser.parse(html, {
                        dt:url_parts.path.split('/')[3],
                        page_url:url,
                        page_hash:htmlMD5
                    });

                    //set new docment
                    var resultNew = new ResultModel(result_json);

                    //save new document
                    crud_service.save(resultNew);
                }
            } else {
                console.log(JSON.stringify(error));
            }
        });
    });
}
module.exports = {
    race:function(meeting){

        console.log('TODO:: parse races/results ' + meeting.track);
        _.filter(meeting.races, function(ele){
            return ele.abandoned !== 'ABD';
        })
        .forEach(function(raceElement, index, array){
            (function() {
                setTimeout(function () {
                    get_runners(conf.sporting_life_com.url + raceElement.racecard_url)
                }, Math.floor(index + Math.random() * conf.sporting_life_com.max_timeout));
            })();
            if (raceElement.result_url !== undefined) {
                (function () {
                    setTimeout(function () {
                        get_result(conf.sporting_life_com.url + raceElement.result_url)
                    }, Math.floor(index + Math.random() * conf.sporting_life_com.max_timeout));
                })();
            }
        });

    }
}