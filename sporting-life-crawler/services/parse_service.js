var _ = require('underscore');

copy_properties_to_parse_obj = function(sourceObj, parseObj){
    _.each(_.keys(sourceObj),function(key){
        parseObj.set(key, sourceObj[key]);
    });
}


var Parse = require('parse/node');
Parse.initialize('test', 'test');
Parse.serverURL = 'http://localhost:1337/parse';


module.exports = {
    save_obj : function(pObjType,pObj,parentCollection,objIndex){
        var obj = new Parse.Object(pObjType);
        copy_properties(pObj, obj);
        if(parentCollection !== null) {
            parentCollection[objIndex] = obj;
        }
        obj.save().then(function (obj) {
            console.log(obj.toJSON());
            // var query = new Parse.Query('Racecard');
            // query.get(obj.id).then(function (objAgain) {
            //     console.log(objAgain.toJSON());
            // }, function (err) {
            //     console.log(err);
            // });
        }, function (err) {
            console.log(err);
        });
    },
    to_parse_obj : function(pObjType,pObj){
        var obj = new Parse.Object(pObjType);
        copy_properties(pObj, obj);
        return obj;
    }
}

