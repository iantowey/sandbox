'use strict';

var neo4j = require('neo4j');

var db = new neo4j.GraphDatabase({
    url: 'http://localhost:7474',
    auth: {
        username: 'neo4j',
        password: 'Lecarr0w'
    }
});

function create_callback(err, results) {
    if (err) throw err;
    console.log("create response  " + JSON.stringify(results) );
};

function query_callback(err, results) {
    if (err) throw err;
    var result = results[0];
    if (!result) {
        console.log('No user found.');
    } else {
        var user = result['user'];
        console.log(user);
    }
};

module.exports = {
    create_tracker:function(properties,post_processing){
        db.cypher({
            query: 'MERGE (p:Person:Tracker { name:{name},account_id:{account_id},dt:{dt} }) RETURN p',
            params: properties,
        }, post_processing);
    },
    create_trainer:function(properties,post_processing){
        db.cypher({
            query: 'MERGE (p:Person:Trainer { name:{name},profile_url:{profile_url},dt:{dt} }) RETURN p',
            params: properties,
        }, post_processing);
    },
    create_jockey:function(properties,post_processing){
        db.cypher({
            query: 'MERGE (p:Person:Jockey { name:{name}, profile_url:{profile_url},dt:{dt} }) RETURN p',
            params: properties,
        }, post_processing);
    },
    create_owner:function(properties,post_processing){
        db.cypher({
            query: 'MERGE (p:Person:Owner { name:{name}, profile_url:{profile_url},dt:{dt} }) RETURN p',
            params: properties,
        }, post_processing);
    },
    create_horse:function(properties,post_processing){
        db.cypher({
            query: 'MERGE (h:Horse:Runner { name:{name}, profile_url:{profile_url},dt:{dt} }) RETURN h',
            params: properties,
        }, post_processing);
    },
    create_sire:function(properties,post_processing){
        db.cypher({
            query: 'MERGE (h:Horse:Sire { name:{name}, profile_url:{profile_url},dt:{dt} }) RETURN h',
            params: properties,
        }, post_processing);
    },
    create_dame:function(properties,post_processing){
        db.cypher({
            query: 'MERGE (h:Horse:Dame { name:{name}, profile_url:{profile_url},dt:{dt} }) RETURN h',
            params: properties,
        }, post_processing);
    },
    create_relationship_horse_dame:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (dame:Horse { name:{dame_name} }),(runner:Horse { name:{runner_name} }) CREATE UNIQUE (dame)-[r:DAME_OF {rel_dt:{rel_dt}}]->(runner)',
            params: properties,
        }, post_processing);
    },
    create_relationship_horse_sire:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (sire:Horse { name:{sire_name} }),(runner:Horse { name:{runner_name} }) CREATE UNIQUE (sire)-[r:SIRE_OF {rel_dt:{rel_dt}}]->(runner)',
            params: properties,
        }, post_processing);
    },
    create_relationship_horse_trainer:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (trainer:Person:Trainer { name:{trainer_name} }),(runner:Horse:Runner { name:{runner_name}}) CREATE (trainer)-[r:TRAINER_OF {rel_dt:{rel_dt}}]->(runner)',
            params: properties,
        }, post_processing);
    },
    create_relationship_horse_owner:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (owner:Person { name:{owner_name} }),(runner:Horse { name:{runner_name} }) CREATE (owner)-[r:TRAINER_OF {rel_dt:{rel_dt}}]->(runner)',
            params: properties,
        }, post_processing);
    },
    create_relationship_horse_jockey:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (jockey:Person:Jockey { name:{jockey_name} }),(runner:Horse:Runner { name:{runner_name} }) CREATE (jockey)-[r:RIDER_OF {rel_dt:{rel_dt}}]->(runner)',
            params: properties,
        }, post_processing);
    },
    create_relationship_tracker_horse:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (tracker:Person { name:{tracker_name} }),(runner:Horse { name:{runner_name} }) CREATE (tracker)-[r:TRACKS {rel_dt:{rel_dt}}]->(runner)',
            params: properties,
        }, post_processing);
    },
    create_relationship_tracker_dame:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (tracker:Person { name:{tracker_name} }),(dame:Horse { name:{dame_name}}) CREATE (tracker)-[r:TRACKS {rel_dt:{rel_dt}}]->(dame)',
            params: properties,
        }, post_processing);
    },
    create_relationship_tracker_sire:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (tracker:Person { name:{tracker_name} }),(sire:Horse { name:{sire_name}}) CREATE (tracker)-[r:TRACKS {rel_dt:{rel_dt}}]->(sire)',
            params: properties,
        }, post_processing);
    },
    create_relationship_tracker_jockey:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (tracker:Person { name:{tracker_name} }),(jockey:Person { name:{jockey_name} }) CREATE (tracker)-[r:TRACKS {rel_dt:{rel_dt}}]->(jockey)',
            params: properties,
        }, post_processing);
    },
    create_relationship_tracker_trainer:function(properties,post_processing){
        db.cypher({
            query: 'MATCH (tracker:Person { name:{tracker_name} }),(trainer:Person { name:{trainer_name}}) CREATE (tracker)-[r:TRACKS {rel_dt:{rel_dt}}]->(trainer)',
            params: properties,
        }, post_processing);
    }
}
