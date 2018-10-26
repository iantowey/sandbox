package com.lucidlogic.horsetracker.model.dto;

import com.lucidlogic.horsetracker.config.Constants;
import com.parse.ParseClassName;
import com.parse.ParseObject;

/**
 * Created by ian on 23/07/16.
 */
    @ParseClassName("Entity")
public class EntityDTO extends ParseObject{

    public String getName(){
        return getString(Constants.NAME);
    }
    public String getProfileUrl(){
        return getString(Constants.PROFILE_URL);
    }
    public String getType(){
        return getString(Constants.TYPE);
    }
    public EntityDetailHistoricalDTO getEntityDetailHistoricalDTO(){ return (EntityDetailHistoricalDTO)get(Constants.ENTITY_DETAIL_HIST); }
    public EntityDetailFutureDTO getEntityDetailFutureDTO(){ return (EntityDetailFutureDTO)get(Constants.ENTITY_DETAIL_FUT); }

    public void setName(String name){ put(Constants.NAME, name) ; }
    public void setProfileUrl(String  profileUrl){ put(Constants.PROFILE_URL, profileUrl); }
    public void setType(String type){
        put(Constants.TYPE, type);
    }
    public void setEntityDetailHistoricalDTO(EntityDetailHistoricalDTO  entityDetailHistoricalDTO){ put(Constants.ENTITY_DETAIL_HIST, entityDetailHistoricalDTO); }
    public void setEntityDetailFutureDTO(EntityDetailFutureDTO entityDetailFutureDTO){ put(Constants.ENTITY_DETAIL_FUT, entityDetailFutureDTO); }


}
