package com.lucidlogic.horsetracker.model;

import android.databinding.BaseObservable;
import android.databinding.Bindable;
import android.os.Parcel;
import android.os.Parcelable;

import com.lucidlogic.horsetracker.BR;

/**
 * Created by ian on 23/07/16.
 */
public class Entity extends BaseObservable implements Parcelable{

    private String id;
    private String type;
    private String name;
    private String profileUrl;
    private boolean userFollowing;
    private EntityDetailFuture entityDetailFuture;
    private EntityDetailHistorical entityDetailHistorical;

    public Entity() {
    }

    protected Entity(Parcel in) {
        id = in.readString();
        type = in.readString();
        name = in.readString();
        profileUrl = in.readString();
        userFollowing = in.readByte() != 0;
        entityDetailFuture = in.readParcelable(this.getClass().getClassLoader());
        entityDetailHistorical = in.readParcelable(this.getClass().getClassLoader());
    }

    @Override
    public void writeToParcel(Parcel dest, int flags) {
        dest.writeString(id);
        dest.writeString(type);
        dest.writeString(name);
        dest.writeString(profileUrl);
        dest.writeByte((byte)(userFollowing ? 1 : 0));
        dest.writeParcelable(entityDetailFuture, 0);
        dest.writeParcelable(entityDetailHistorical, 0);

    }

    @Override
    public int describeContents() {
        return 0;
    }

    public static final Creator<Entity> CREATOR = new Creator<Entity>() {
        @Override
        public Entity createFromParcel(Parcel in) {
            return new Entity(in);
        }

        @Override
        public Entity[] newArray(int size) {
            return new Entity[size];
        }
    };

    @Bindable
    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
        notifyPropertyChanged(BR.id);
    }

    @Bindable
    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
        notifyPropertyChanged(BR.type);

    }

    @Bindable
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
        notifyPropertyChanged(BR.name);

    }

    @Bindable
    public String getProfileUrl() {
        return profileUrl;

    }

    public void setProfileUrl(String profileUrl) {
        this.profileUrl = profileUrl;
        notifyPropertyChanged(BR.profileUrl);
    }

    @Bindable
    public EntityDetailFuture getEntityDetailFuture() {
        return entityDetailFuture;
    }

    public void setEntityDetailFuture(EntityDetailFuture entityDetailFuture) {
        this.entityDetailFuture = entityDetailFuture;
        notifyPropertyChanged(BR.entityDetailFuture);

    }

    @Bindable
    public EntityDetailHistorical getEntityDetailHistorical() {
        return entityDetailHistorical;
    }

    public void setEntityDetailHistorical(EntityDetailHistorical entityDetailHistorical) {
        this.entityDetailHistorical = entityDetailHistorical;
        notifyPropertyChanged(BR.entityDetailHistorical);

    }

    @Bindable
    public boolean getUserFollowing() {
        return userFollowing;
    }

    public void setUserFollowing(boolean userFollowing) {
        this.userFollowing = userFollowing;
        notifyPropertyChanged(BR.userFollowing);

    }
}
