var express = require('express');
var router = express.Router();
var ResultModel = require("../../models/result");

router.get('/', function(req, res, next) {
    ResultModel.find({}, function(error, results) {
        if(error) {
            console.log("An error happened -> " + JSON.stringify(error));
        } else {
            res.send(JSON.stringify(results));
        }
    });
});

router.get("/id", function(req, res, next) {
    var pageUrl =  req.query.q;
    console.log("result url  -> " + pageUrl);
    ResultModel.findByPageUrl(pageUrl, function(error, results) {
        if(error) {
            console.log("An error happened -> " + JSON.stringify(error));
        } else {
            res.send(JSON.stringify(results));
        }
    });
});

module.exports = router;
