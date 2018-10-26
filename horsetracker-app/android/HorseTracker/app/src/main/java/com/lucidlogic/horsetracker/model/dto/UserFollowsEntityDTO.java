package com.lucidlogic.horsetracker.model.dto;

import com.lucidlogic.horsetracker.config.Constants;
import com.parse.ParseClassName;
import com.parse.ParseObject;
import com.parse.ParseUser;

/**
 * Created by itowey on 12/08/16.
 */
@ParseClassName(value = "UserFollowsEntity")
public class UserFollowsEntityDTO extends ParseObject{

    public ParseUser getUser(){return getParseUser(Constants.USER);}
    public EntityDTO getEntity(){return (EntityDTO)getParseObject(Constants.ENTITY);}

    public void setUser(ParseUser user){put(Constants.USER, user);}
    public void setEntity(EntityDTO entity){put(Constants.ENTITY, entity);}

}
