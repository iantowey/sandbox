var express = require('express');
var router = express.Router();
var RacecardModel = require("../../models/racecard");

/* GET racecard list . */
router.get('/', function(req, res, next) {
  RacecardModel.find({}, function(error, results) {
    if(error) {
      console.log("An error happened -> " + JSON.stringify(error));
    } else {
      res.send(JSON.stringify(results));
    }
  });
});

router.get(/^\/([0-9]{2}-[0-9]{2}-[0-9]{4})$/, function(req, res, next) {
  var dt = req.params[0];
  RacecardModel.findByDate(dt, function(error, results) {
    if(error) {
      console.log("An error happened -> " + JSON.stringify(error));
    } else {
      res.send(JSON.stringify(results));
    }
  });
});

router.get('/crawl', function(req, res, next) {
  var dt = req.query.dt;
  RacecardModel.findByDate(dt, function(error, results) {
    if(error) {
      console.log("An error happened -> " + JSON.stringify(error));
    } else {
      res.send(JSON.stringify(results));
    }
  });
});

module.exports = router;
