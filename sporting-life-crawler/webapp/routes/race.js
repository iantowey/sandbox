var express = require('express');
var router = express.Router();
var RaceModel = require("../../models/race");

router.get('/', function(req, res, next) {
    RaceModel.find({}, function(error, results) {
        if(error) {
            console.log("An error happened -> " + JSON.stringify(error));
        } else {
            res.send(JSON.stringify(results));
        }
    });
});

router.get("/id", function(req, res, next) {
    var pageUrl =  req.query.q;
    console.log("race url  -> " + pageUrl);
    RaceModel.findByPageUrl(pageUrl, function(error, results) {
        if(error) {
            console.log("An error happened -> " + JSON.stringify(error));
        } else {
            res.send(JSON.stringify(results));
        }
    });
});

module.exports = router;
