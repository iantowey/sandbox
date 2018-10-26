package com.lucidlogic.horsetracker.handler;

import com.annimon.stream.Stream;
import com.lucidlogic.horsetracker.config.Constants;
import com.lucidlogic.horsetracker.model.Entity;
import com.lucidlogic.horsetracker.model.dto.EntityDTO;
import com.lucidlogic.horsetracker.model.dto.UserFollowsEntityDTO;
import com.lucidlogic.horsetracker.service.ParseService;
import com.parse.DeleteCallback;
import com.parse.ParseException;
import com.parse.ParseObject;
import com.parse.ParseQuery;
import com.parse.ParseUser;
import com.parse.SaveCallback;

import java.util.List;

import timber.log.Timber;

/**
 * Created by ian on 02/08/16.
 */
public class EntityClickHandler {

    public void onEntityClick(Entity e){
        Timber.i("%s %s clicked" ,e.getName(), e.getId());
        if(!e.getUserFollowing()) {
            e.setUserFollowing(true);
            UserFollowsEntityDTO userFollowsEntityDTO = ParseObject.create(UserFollowsEntityDTO.class);
            userFollowsEntityDTO.setUser(ParseUser.getCurrentUser());
            EntityDTO entityDTO = ParseObject.createWithoutData(EntityDTO.class, e.getId());
            userFollowsEntityDTO.setEntity(entityDTO);
            userFollowsEntityDTO.saveInBackground(new SaveCallback() {
                @Override
                public void done(ParseException e) {

                    if (e == null) {
                        Timber.i("userFollowsEntityDTO saved");
                    } else {
                        e.printStackTrace();
                    }
                }
            });
        } else {
            e.setUserFollowing(false);
            try {
                EntityDTO entityDTO = ParseObject.createWithoutData(EntityDTO.class, e.getId());
                List<UserFollowsEntityDTO> userFollowsEntityDTOs =
                        ParseQuery.getQuery(UserFollowsEntityDTO.class)
                                .whereEqualTo(Constants.USER, ParseUser.getCurrentUser())
                                .whereEqualTo(Constants.ENTITY, entityDTO).find();
                Stream.of(userFollowsEntityDTOs).forEach(userFollowsEntityDTO -> userFollowsEntityDTO.deleteEventually(new DeleteCallback() {
                    @Override
                    public void done(ParseException e) {
                        Timber.i("Removed tracking");
                    }
                }));
                ParseService.getEntitiesUserTracks();
            } catch (ParseException e1) {
                e1.printStackTrace();
            }
        }
    }
}
