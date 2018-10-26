package com.lucidlogic.horsetracker.service;

import com.lucidlogic.horsetracker.config.Constants;
import com.lucidlogic.horsetracker.model.dto.EntityDTO;
import com.lucidlogic.horsetracker.model.dto.RaceDTO;
import com.lucidlogic.horsetracker.model.dto.RacecardDTO;
import com.lucidlogic.horsetracker.model.dto.UserFollowsEntityDTO;
import com.parse.FindCallback;
import com.parse.ParseException;
import com.parse.ParseObject;
import com.parse.ParseQuery;
import com.parse.ParseUser;

import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.List;

import rx.Observable;
import rx.parse.ParseObservable;

/**
 * Created by ian on 10/07/16.
 */
public class ParseService {

    public static Observable<RacecardDTO> getRacecards(Date raceDate){

        Calendar cal = Calendar.getInstance(); // locale-specific
        cal.setTime(raceDate);
        cal.set(Calendar.HOUR_OF_DAY, 0);
        cal.set(Calendar.MINUTE, 0);
        cal.set(Calendar.SECOND, 0);
        cal.set(Calendar.MILLISECOND, 0);
        long time = cal.getTimeInMillis();
        return ParseObservable.find(ParseQuery.getQuery(RacecardDTO.class)
                .include("meetings")
//                .whereEqualTo(Constants.DATE, cal.getTime())
                .include("meetings.races")
                .addAscendingOrder("date"))
                .first();
    }

    public static Observable<RaceDTO> getRace(String raceObjectId) {

        return ParseObservable.find(
                ParseQuery.getQuery(RaceDTO.class)
                .whereEqualTo(Constants.OBJECT_ID, raceObjectId)
                .include(Constants.RUNNERS)
                .include(Constants.RUNNERS + "." + Constants.JOCKEY)
                .include(Constants.RUNNERS + "." + Constants.HORSE)
                .include(Constants.RUNNERS + "." + Constants.TRAINER)
        );
    }

    public static <T extends ParseObject> void saveToLocalDb(String label, T parseObject) throws ParseException {
        List<T> parseObjects = new ArrayList<>();
        parseObjects.add(parseObject);
        saveToLocalDb(label, parseObjects);
    }

    public static <T extends ParseObject> void saveToLocalDb(String label, List<T> parseObjects) throws ParseException {
        ParseObject.pinAll(label, parseObjects);
    }

    public static void removeFreomLocalDb(String label) throws ParseException {
        ParseObject.unpinAll(label);

    }

    public static Observable<EntityDTO> getEntity(String entityObjectId) {

        return ParseObservable.find(
                ParseQuery.getQuery(EntityDTO.class)
                        .whereEqualTo(Constants.OBJECT_ID, entityObjectId)
                        .include(Constants.ENTITY_DETAIL_FUT)
                        .include(Constants.ENTITY_DETAIL_HIST));
    }

    public static void putInLocalParseStore(ParseObject parseObject){
        parseObject.pinInBackground();
    }

    public static <T extends ParseObject> T getFromLocalParseStore(String objectId, Class<T> clazz){
        try {
            return ParseQuery.getQuery(clazz).fromLocalDatastore().get(objectId);
        } catch (ParseException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static <T extends ParseObject> List<T> getFromLocalParseStoreByLabel(String objectLabel, Class<T> clazz){
        try {
            return ParseQuery.getQuery(clazz).fromLocalDatastore().fromPin(objectLabel).find();
        } catch (ParseException e) {
            e.printStackTrace();
            return null;
        }
    }

    public static void initUser() {
        getEntitiesUserTracks();
    }


    public static void getEntitiesUserTracks() {
        ParseUser parseUser = ParseUser.getCurrentUser();
        if(parseUser != null){
            //TODO: get all user saved data
            ParseQuery<UserFollowsEntityDTO> trackEntityDTOParseQuery = ParseQuery.getQuery(UserFollowsEntityDTO.class);
            trackEntityDTOParseQuery.whereEqualTo("user", parseUser).findInBackground(new FindCallback<UserFollowsEntityDTO>() {
                @Override
                public void done(List<UserFollowsEntityDTO> objects, ParseException e) {
                    try {
                        saveToLocalDb(Constants.USER_FOLLOWS_ENTITY_LIST, objects);
                    } catch (ParseException e1) {
                        e1.printStackTrace();
                    }
                }
            });

        }
    }
}
