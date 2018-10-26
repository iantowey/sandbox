var Parse = require('parse/node');
Parse.initialize('test', 'test');
Parse.serverURL = 'http://localhost:1337/parse';

var subObj = new Parse.Object('GameScore');


var obj = new Parse.Object('GameScore');
obj.set('score', 1337);
obj.save().then(function (obj) {
    console.log(obj.toJSON());
    var query = new Parse.Query('GameScore');
    query.get(obj.id).then(function (objAgain) {
        console.log(objAgain.toJSON());
    }, function (err) {
        console.log(err);
    });
}, function (err) {
    console.log(err);
});
