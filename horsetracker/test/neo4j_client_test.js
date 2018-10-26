var neo4j_client = require("../services/integrations/neo4j_client");
var chai   = require("chai");
var fs = require('fs');
var file_name = './tmp/detail_json_test';

describe("Parse horse details", function() {


    before(function(done){
        setTimeout(function(){
            foo = true;
            done();
        }, 1000);
    });

    var horse_profile = JSON.parse(fs.readFileSync(file_name));
    console.log(JSON.stringify(horse_profile));

    it("Create horse Node test", function (done) {

        neo4j_client.create_horse({
            name:horse_profile.name.trim(),
            age:horse_profile.age,
            sex:horse_profile.sex,
            profile_url:"TODO"
        },done);
    });

    it("Create Sire Node test", function (done) {

        neo4j_client.create_sire({
            name:horse_profile.sire.name.trim(),
            profile_url:horse_profile.sire.profile_url
        },done);
    });

    it("Create Dame Node test", function (done) {

        neo4j_client.create_dame({
            name:horse_profile.dam.name.trim(),
            profile_url:horse_profile.dam.profile_url
        },done);
    });

    it("Create trainer Node test", function (done) {

        neo4j_client.create_trainer({
            name:horse_profile.trainer.name.trim(),
            profile_url:horse_profile.trainer.profile_url
        },done);
    });

    it("Create owner Node test", function (done) {
v 
        neo4j_client.create_owner({
            name:horse_profile.owner.trim(),
            profile_url:""
        },done);
    });

    it("Create tracker Node test", function (done) {

        neo4j_client.create_tracker({
            name:'Ian Towey',
            account_id:"16234"
        },done);
    });

    it("Create runner/sire relationship", function (done) {

        neo4j_client.create_relationship_horse_sire({
            sire_name:horse_profile.sire.name.trim(),
            runner_name:horse_profile.name
        },done);
    });

    it("Create runner/dame relationship", function (done) {

        neo4j_client.create_relationship_horse_dame({
            dame_name:horse_profile.dam.name.trim(),
            runner_name:horse_profile.name
        },done);
    });

    it("Create runner/owner relationship", function (done) {

        neo4j_client.create_relationship_horse_owner({
            owner_name:horse_profile.owner.trim(),
            runner_name:horse_profile.name
        },done);
    });

    it("Create runner/trainer relationship", function (done) {

        neo4j_client.create_relationship_horse_trainer({
            trainer_name:horse_profile.trainer.name.trim(),
            runner_name:horse_profile.name
        },done);
    });

    it("Create tracker/horse relationship", function (done) {

        neo4j_client.create_relationship_tracker_horse({
            tracker_name:horse_profile.sire.name.trim(),
            runner_name:'Ian Towey'
        },done);
    });

});