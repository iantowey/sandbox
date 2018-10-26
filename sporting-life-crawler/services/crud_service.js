
module.exports = {
    save:function(obj){
        obj.save(function(error) {
            if(error) {
                console.error(error);
            } else {
                console.log(obj.constructor.name + " saved");
            }
        });
    },
    update: function(obj){

    },
    delete:function(obj){
        obj.remove(function (error) {
            if (error) {
                console.log("An error happened -> " + JSON.stringify(error));
            }
            console.log(obj.constructor.name  + " document was removed!");
        });
    },
    select:function(obj){
        throw new Error;
    }
}